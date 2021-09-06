//
//  MatrixPivotFinder.swift
//  SwiftyMath
//
//  Created by Taketo Sano on 2019/10/26.
//

//  Implementation based on:
//
//  "Parallel Sparse PLUQ Factorization modulo p", Charles Bouillaguet, Claire Delaplace, Marie-Emilie Voge.
//  https://hal.inria.fr/hal-01646133/document
//
//  see also:
//
//  SpaSM (Sparse direct Solver Modulo p)
//  https://github.com/cbouilla/spasm

import SwmCore
import SwmMatrixTools
import Dispatch

public final class _MatrixPivotFinder {
    public let mode: PivotMode
    public let size: MatrixSize
    
    private let data: [[Int]]
    private let rowHeads: [Int: Int]
    private let rowWeights: [Double]
    private let colWeights: [Double]
    private let candidates: [Int: Set<Int>]
    
    private var pivotRows: Set<Int> // Set(rows)
    private var pivotTable: [Int : Int] // [col : row]
    private var result: [(Int, Int)]

    public var debug: Bool = false
    
    public init<Impl: MatrixImpl>(_ A: Impl, mode: PivotMode = .rowBased, pivotWeightThreshold: Double = 1.0) where Impl.BaseRing: ComputationalRing {
        let (n, m) = (mode == .rowBased)
            ? A.size
            : (A.size.cols, A.size.rows)
        
        var data: [[Int]] = Array(repeating: [], count: n)
        var rowHeads: [Int: Int] = [:]
        var rowWeights: [Double] = Array(repeating: 0, count: n)
        var colWeights: [Double] = Array(repeating: 0, count: m)
        var candidates: [Int: Set<Int>] = [:]

        for (i0, j0, a) in A.nonZeroEntries {
            let (i, j) = (mode == .rowBased) ? (i0, j0) : (j0, i0)
            let w = a.computationalWeight
            
            data[i].append(j)
            rowWeights[i] += w
            colWeights[j] += w
            
            if rowHeads[i] == nil || j < rowHeads[i]! {
                rowHeads[i] = j
            }
            if a.isInvertible && a.computationalWeight <= pivotWeightThreshold {
                candidates[i, default: []].insert(j)
            }
        }
        
        self.mode = mode
        self.size = A.size
        self.data = data
        self.rowHeads = rowHeads
        self.rowWeights = rowWeights
        self.colWeights = colWeights
        self.candidates = candidates
        self.pivotTable = [:]
        self.pivotRows = []
        self.result = []
    }
    
    public convenience init<Impl, n, m>(_ A: MatrixIF<Impl, n, m>, mode: PivotMode = .rowBased) where Impl.BaseRing: ComputationalRing {
        self.init(A.impl, mode: mode)
    }
    
    public var pivots: [(Int, Int)] {
        mode == .rowBased
            ? result
            : result.map{ (i, j) in (j, i) }
    }
    
    public var rowPermutation: Permutation<anySize> {
        asPermutation(size.rows, pivots.map{ $0.0 })
    }
    
    public var colPermutation: Permutation<anySize> {
        asPermutation(size.cols, pivots.map{ $0.1 })
    }
    
    private func asPermutation(_ length: Int, _ order: [Int]) -> Permutation<anySize> {
        Permutation.fill(length: length, indices: order).inverse!
    }
    
    public func run() {
        log("Start pivot search.")
        log("size: \(size)")
        
        findFLPivots()
        findFLColumnPivots()
        findCycleFreePivots()
        
        let pivots = sortPivots()
        
        if pivots.isEmpty {
            log("No pivots found.")
        } else {
            self.result = pivots
            
            log("Total: \(pivots.count)\(pivots.allSatisfy({ (i, _) in data[i].count == 1 }) ? " perfect" : "") pivots.")
        }
    }
    
    // Faugère-Lachartre pivot search
    private func findFLPivots() {
        var candidates: [Int : Int] = [:] // col -> row
        
        for (i, j) in rowHeads where isCandidate(i, j) {
            if !candidates.contains(key: j) || weight(i, j) < weight(candidates[j]!, j) {
                candidates[j] = i
            }
        }
        
        candidates.forEach{ (j, i) in
            setPivot(i, j)
        }
        
        log("FL-pivots: \(candidates.count)")
    }
    
    private func findFLColumnPivots() {
        if pivotRows.count >= min(size.rows, size.cols) {
            return
        }
        
        let count = pivotTable.count
        var occupiedCols = Set(pivotRows.flatMap{ i in data[i] })
        
        for i in data.indices where !pivotRows.contains(i) {
            var currentCols: Set<Int> = []
            var candidates: [Int] = []
            
            for j in data[i] where !occupiedCols.contains(j) {
                currentCols.insert(j)
                if isCandidate(i, j) {
                    candidates.append(j)
                }
            }
            
            if let j = candidates.min(by: { j in weight(i, j) }) {
                setPivot(i, j)
                occupiedCols.formUnion(currentCols)
            }
        }
        
        log("FL-col-pivots: \(pivotTable.count - count)")
    }
    
    private func findCycleFreePivots() {
        var remainingRows = Set(data.indices).subtracting(pivotRows).sorted{ -rowWeights[$0] } // sorted in heavy order.
        if remainingRows.isEmpty {
            return
        }
        
//        log("start findCycleFreePivots.")
        
        let count = pivotTable.count
        var batch = 1
        var itr = 1
        
        while !remainingRows.isEmpty {
//            log("iteration: \(itr), targets: \(min(batch, remainingRows.count)) / \(remainingRows.count)")
            
            let targets = (0 ..< batch).compactMap { _ in
                remainingRows.isEmpty ? nil : remainingRows.removeLast()
            }
            let candidates = targets.parallelCompactMap { i in
                self.findCycleFreePivot(inRow: i) // runs in parallel
            }
            
//            log("found \(candidates.count) candidates")
            
            if candidates.count == 0 {
                batch *= 2
            } else if candidates.count >= 2 && batch > 1 {
                batch /= 2
            }
            
            if let (i, j) = candidates.min(by: { (i, j) in weight(i, j) }) {
                self.setPivot(i, j)
                for (i1, _) in candidates where i1 != i {
                    remainingRows.append(i1)
                }
            }
            itr += 1
        }

        log("cycle-free-pivots: \(pivotTable.count - count)")
    }

    private func findCycleFreePivot(inRow i: Int) -> (Int, Int)? {
        let row = data[i]

        //       j1
        //  i [  O     O     C  C     ]    O: queued, C: candidates
        //       |     |     ^
        //       V     |     | rmv
        // i2 [  X     o     o     O  ]
        //             |           |
        //             V           |
        // i3 [        X           |  ]
        //                         |

        // the following are col-indexed.
        var candidates: Set<Int> = []
        var queue: OrderedSet<Int> = []

        // initialize
        for j in row {
            if pivotTable.contains(key: j) {
                queue.append(j)
            } else if isCandidate(i, j) {
                candidates.insert(j)
            }
        }
        
        var idx = 0
        while idx < queue.count && !candidates.isEmpty {
            let j1 = queue[idx]
            let i2 = pivotTable[j1]!

            for j2 in data[i2] {
                if pivotTable.contains(key: j2) {
                    queue.append(j2) // will be ignored if already queued.
                } else if candidates.contains(j2) {
                    candidates.remove(j2)
                    if candidates.isEmpty {
                        break
                    }
                }
            }
            
            idx += 1
        }

        if let j = candidates.min(by: { j in weight(i, j) }) {
            return (i, j)
        } else {
            return nil
        }
    }
    
    private func sortPivots() -> [(Int, Int)] {
        try! pivotTable.keys.topologicalSort { j1 in
            let i1 = pivotTable[j1]!
            return data[i1].filter { j2 in
                j1 != j2 && pivotTable.contains(key: j2)
            }
        }.map { j in
            (pivotTable[j]!, j)
        }
    }
    
    private func isCandidate(_ i: Int, _ j: Int) -> Bool {
        candidates[i]?.contains(j) ?? false
    }
    
    private func weight(_ i: Int, _ j: Int) -> Double {
        rowWeights[i] * colWeights[j]
    }
    
    private func setPivot(_ i: Int, _ j: Int) {
        pivotRows.insert(i)
        pivotTable[j] = i
    }
    
    private func log(_ msg: @autoclosure () -> String) {
        if debug {
            print(msg())
        }
    }
}
