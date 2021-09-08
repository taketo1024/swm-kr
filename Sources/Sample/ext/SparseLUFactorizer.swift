//
//  LUFactorizer.swift
//  SwiftyMath
//
//  Created by Taketo Sano on 2019/10/29.
//

import Foundation
import SwmCore
import SwmMatrixTools

public final class _SparseLUFactorizer<M: SparseMatrixImpl & LUFactorizable> where M.BaseRing: ComputationalRing {
    public typealias Matrix = M
    
    public private(set) var target: Matrix
    public private(set) var P: Permutation<anySize>
    public private(set) var Q: Permutation<anySize>
    public private(set) var L: Matrix
    public private(set) var U: Matrix
    
    public var densityThreshold = 0.3
    public var debug: Bool
    
    public convenience init<n, m>(_ A: MatrixIF<Matrix, n, m>, debug: Bool = false) {
        self.init(A.impl, debug: debug)
    }
    
    public init(_ A: Matrix, debug: Bool = false) {
        let (n, m) = A.size
        self.target = A
        self.P = Permutation.identity(length: n)
        self.Q = Permutation.identity(length: m)
        self.L = Matrix.zero(size: (n, 0))
        self.U = Matrix.zero(size: (0, m))
        self.debug = debug
    }
    
    public var isComplete: Bool {
        target.isZero
    }
    
    public var result: (P: Permutation<anySize>, Q: Permutation<anySize>, L: Matrix, U: Matrix) {
        (P, Q, L, U)
    }
    
    public func run() {
        if target.isZero {
            return
        }
        
        log("start SparseLUFactorization, target: \(target.size), density: \(target.density)")

        var itr = 1
        var pivotWeight = 1.0
        
        while !target.isZero && target.density < densityThreshold {
            log("partialLU: \(itr), target: \(target.size), density: \(target.density), pivotWeight: \(pivotWeight)")
            
            let r = partialLU()
            if r == 0 {
                if pivotWeight == 1.0 {
                    pivotWeight = Double.infinity
                } else {
                    break
                }
            }
            
            itr += 1
        }
        
        if !target.isZero {
            log("fullLU, target: \(target.size), density: \(target.density)")
            fullLU()
        }
        
        log("Done, rank: \(L.size.cols)")
        log("")
    }
    
    @discardableResult
    public func partialLU() -> Int {
        // If
        //
        //   PAQ = [L, B]
        //         [C, D]
        //
        // with L: lower-triangular, then
        //
        //   PAQ = [L, 0] * [I, B']
        //         [C, S]   [0, I ]
        //
        //       = [L] * [I, B'] + [0, 0]
        //         [C]             [0, S] .
        //
        // where
        //
        //    B' = L^-1 B,
        //    S  = D - C L^-1 B.
        //

        let (n, m) = target.size
        let isLarge = min(n, m) > 20000
        
        let pf = _MatrixPivotFinder(target, mode: .colBased)
        pf.run()
        
        let r = pf.pivots.count
        
        if r == 0 {
            log("no pivots found.")
            return 0
        } else {
            log("found \(r) pivots")
        }

        let (P, Q) = (pf.rowPermutation, pf.colPermutation)

        let R  = target.permute(rowsBy: pf.rowPermutation, colsBy: pf.colPermutation)
        let L0 = R.submatrix(rowRange: 0 ..< r, colRange: 0 ..< r)
        let B  = R.submatrix(rowRange: 0 ..< r, colRange: r ..< m)
        let C  = R.submatrix(rowRange: r ..< n, colRange: 0 ..< r)
        let D  = R.submatrix(rowRange: r ..< n, colRange: r ..< m)
        
        if isLarge {
            log("solve lower triangular")
        }
        
        let B1 = Matrix.solveLowerTriangular(L0, B)
        
        let L = L0.stack(C)
        let U = Matrix.identity(size: (r, r)).concat(B1)
        
        compose(P, Q, L, U)
        
        if isLarge {
            log("compute S")
        }
        
        let S = D - C * B1
        self.target = S
        
        return r
    }
    
    public func fullLU() {
        if target.isZero {
            return
        }
        
        let f = LUFactorizer(target)
        f.run()
        let (P, Q, L, U) = f.result
        
        compose(P, Q, L, U)
        target = .zero(size: (L.size.rows - L.size.cols, U.size.cols - U.size.rows))
    }
    
    private func compose(_ P2s: Permutation<anySize>, _ Q2s: Permutation<anySize>, _ L2s: Matrix, _ U2s: Matrix) {
        if L.size.cols == 0 {
            (P, Q, L, U) = (P2s, Q2s, L2s, U2s)
            return
        }
        
        let r1 = L.size.cols
        let r2 = L2s.size.cols
        
        let P2 = P2s.shifted(r1)
        let Q2 = Q2s.shifted(r1)
        let L2 = Matrix.zero(size: (r1, r2)).stack (L2s)
        let U2 = Matrix.zero(size: (r2, r1)).concat(U2s)
        
        P = P2 * P
        Q = Q2 * Q
        L = L.permuteRows(by: P2).concat(L2)
        U = U.permuteCols(by: Q2).stack (U2)
    }
    
    private func log(_ msg: @autoclosure () -> String) {
        if debug {
            let df = DateFormatter()
            df.timeStyle = .medium
            df.dateStyle = .short
            df.locale = Locale(identifier: "ja_JP")
            let now = df.string(from: Date())
            print("[\(now)]", msg())
        }
    }
}
