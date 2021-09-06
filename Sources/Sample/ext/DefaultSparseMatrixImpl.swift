//
//  CSCMatrixImpl.swift
//  
//
//  Created by Taketo Sano on 2021/06/04.
//

import SwmCore
import SwmMatrixTools

//  CSC (compressed sparse colums) format.
//  https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)

public struct _DefaultSparseMatrixImpl<R: Ring>: SparseMatrixImpl {
    public typealias BaseRing = R
    
    public let size: (rows: Int, cols: Int)
    
    private let values: [R]        // length: nnz
    private let rowIndices: [Int]  // length: nnz
    private let indexRanges: [Int] // length: m+1
    
    private init(size: MatrixSize, values: [R], rowIndices: [Int], indexRanges: [Int]) {
        assert(!values.contains{ $0.isZero })
        assert(values.count == rowIndices.count)
        assert(indexRanges.count == size.cols + 1 && indexRanges[0] == 0 && indexRanges.last! == values.count)

        self.size = size
        self.values = values
        self.rowIndices = rowIndices
        self.indexRanges = indexRanges
    }
    
    fileprivate init(size: MatrixSize, compressing cols: [[ColEntry<R>]], sorted: Bool) {
        let (values, rowIndices, indexRanges) = Self.compress(cols, sorted: sorted)
        self.init(
            size: size,
            values: values,
            rowIndices: rowIndices,
            indexRanges: indexRanges
        )
    }
    
    @_specialize(where R == 𝐙)
    @_specialize(where R == 𝐐)
    @_specialize(where R == 𝐅₂)
    
    public init(size: MatrixSize, initializer: (Initializer) -> Void) {
        assert(size.rows >= 0)
        assert(size.cols >= 0)
        
        let (n, m) = size
        var cols: [[ColEntry<R>]] = .init(repeating: [], count: m)
        
        initializer { (i, j, a) in
            assert( 0 <= i && i < n )
            assert( 0 <= j && j < m )
            if !a.isZero {
                cols[j].append((i, a))
            }
        }
        
        self.init(size: size, compressing: cols, sorted: false)
    }

    public subscript(i: Int, j: Int) -> R {
        get {
            let r = indexRange(j)
            if let idx = rowIndices[r].firstIndex(of: i) {
                return values[idx]
            } else {
                return .zero
            }
        } set {
            fatalError("not supported")
        }
    }
    
    private func indexRange(_ j: Int) -> Range<Int> {
        indexRanges[j] ..< indexRanges[j + 1]
    }
    
    public var numberOfNonZeros: Int {
        values.count
    }
    
    public var nonZeroEntries: AnySequence<MatrixEntry<R>> {
        AnySequence(NonZeroEntryIterator(self))
    }
    
    public func colVector(_ j: Int) -> Self {
        let r = indexRange(j)
        return .init(
            size: (size.rows, 1),
            values: Array(values[r]),
            rowIndices: Array(rowIndices[r]),
            indexRanges: [0, r.upperBound - r.lowerBound]
        )
    }
    
    public func firstEntry(inCol j: Int = 0) -> ColEntry<R>? {
        let r = indexRange(j)
        if !r.isEmpty {
            let idx = r.lowerBound
            return (rowIndices[idx], values[idx])
        } else {
            return nil
        }
    }
    
    public func lastEntry(inCol _j: Int? = nil) -> ColEntry<R>? {
        let j = _j ?? size.cols - 1
        let r = indexRange(j)
        if !r.isEmpty {
            let idx = r.upperBound - 1
            return (rowIndices[idx], values[idx])
        } else {
            return nil
        }
    }
    
    public static func ==(a: Self, b: Self) -> Bool {
        (a.values, a.rowIndices, a.indexRanges) == (b.values, b.rowIndices, b.indexRanges)
    }
    
    @_specialize(where R == 𝐙)
    @_specialize(where R == 𝐐)
    @_specialize(where R == 𝐅₂)
    
    public static func +(a: Self, b: Self) -> Self {
        assert(a.size == b.size)
        
        let m = a.size.cols
        
        var values: [R] = []
        var rowIndices: [Int] = []
        var indexRanges:  [Int] = []
        
        values.reserveCapacity(a.numberOfNonZeros + b.numberOfNonZeros)
        rowIndices.reserveCapacity(a.numberOfNonZeros + b.numberOfNonZeros)
        indexRanges.reserveCapacity(m + 1)

        var c = 0
        indexRanges.append(0)
        
        for j in 0 ..< m {
            let r1 = a.indexRange(j)
            let r2 = b.indexRange(j)
            
            var idx1 = r1.lowerBound
            var idx2 = r2.lowerBound
            
            while idx1 < r1.upperBound && idx2 < r2.upperBound {
                let i1 = a.rowIndices[idx1]
                let i2 = b.rowIndices[idx2]
                
                if i1 < i2 {
                    rowIndices.append(i1)
                    values.append(a.values[idx1])
                    idx1 += 1
                    c += 1
                } else if i1 > i2 {
                    rowIndices.append(i2)
                    values.append(b.values[idx2])
                    idx2 += 1
                    c += 1
                } else {
                    let x = a.values[idx1] + b.values[idx2]
                    if !x.isZero {
                        rowIndices.append(i1)
                        values.append(x)
                        c += 1
                    }
                    idx1 += 1
                    idx2 += 1
                }
            }
            
            while idx1 < r1.upperBound {
                rowIndices.append(a.rowIndices[idx1])
                values.append(a.values[idx1])
                idx1 += 1
                c += 1
            }
            
            while idx2 < r2.upperBound {
                rowIndices.append(b.rowIndices[idx2])
                values.append(b.values[idx2])
                idx2 += 1
                c += 1
            }
            
            indexRanges.append(c)
        }
        
        return .init(
            size: a.size,
            values: values,
            rowIndices: rowIndices,
            indexRanges: indexRanges
        )
    }
    
    @_specialize(where R == 𝐙)
    @_specialize(where R == 𝐐)
    @_specialize(where R == 𝐅₂)
    
    public static func *(a: Self, b: Self) -> Self {
        assert(a.size.cols == b.size.rows)
        
        //      k              j        j
        // i |  a    *  |     | |    i |*|
        //   |          |   k |b|      | |
        //   |  *       |  x  | |  ->  | |
        //   |          |     |*|      |*|
        //   |  *       |     | |      | |
        //
        
        let cols: [[ColEntry<R>]] =
            Array(0 ..< b.size.cols).parallelMap { j -> [ColEntry<R>] in
                b.indexRange(j).flatMap { idx2 -> [ColEntry<R>] in
                    let k = b.rowIndices[idx2]
                    let b_kj = b.values[idx2]
                    return a.indexRange(k).map { idx1 in
                        let i = a.rowIndices[idx1]
                        let a_ik = a.values[idx1]
                        return ColEntry(i, a_ik * b_kj)
                    }
                }
                .group{ $0.row }
                .mapValues { $0.sum{ $0.value } }
                .compactMap{ (i, a) in a.isZero ? nil : ColEntry(i, a) }
                .sorted{ $0.row }
            }
        
        return .init(
            size: (a.size.rows, b.size.cols),
            compressing: cols,
            sorted: true
        )
    }
    
    public func serialize() -> [R] {
        let (n, m) = size
        var res: [R] = .init(repeating: .zero, count: n * m)
        for (i, j, a) in nonZeroEntries {
            res[i * m + j] = a
        }
        return res
    }
    
    @_specialize(where R == 𝐙)
    @_specialize(where R == 𝐐)
    @_specialize(where R == 𝐅₂)
    
    private static func compress(_ cols: [[ColEntry<R>]], sorted: Bool) -> ([R], [Int], [Int]) {
        let nnz = cols.sum{ $0.count }
        
        var values: [R] = []
        var rowIndices: [Int] = []
        var indexRanges:  [Int] = []
        
        values.reserveCapacity(nnz)
        rowIndices.reserveCapacity(nnz)
        indexRanges.reserveCapacity(cols.count + 1)

        var r = 0
        indexRanges.append(0)
        
        for col in cols {
            let scol = sorted ? col : col.sorted(by: { $0.row })
            for (i, a) in scol {
                values.append(a)
                rowIndices.append(i)
            }
            r += col.count
            indexRanges.append(r)
        }
        
        return (values, rowIndices, indexRanges)
    }

    public struct NonZeroEntryIterator: Sequence, IteratorProtocol {
        public typealias Element = MatrixEntry<R>
        
        private let data: _DefaultSparseMatrixImpl<R>
        private var idx: Int
        private var col: Int
        private let maxCol: Int
        
        public init(_ data: _DefaultSparseMatrixImpl<R>) {
            let m = data.size.cols
            self.data = data
            self.idx = 0
            self.col = data.indexRanges.firstIndex{ $0 > 0 }.flatMap{ $0 - 1 } ?? m
            self.maxCol = m
        }
        
        public mutating func next() -> MatrixEntry<R>? {
            guard idx < data.values.count else {
                return nil
            }
            
            defer {
                idx += 1
                while col < maxCol && idx >= data.indexRanges[col + 1] {
                    col += 1
                }
            }
            
            return (data.rowIndices[idx], col, data.values[idx])
        }
    }
}

extension _DefaultSparseMatrixImpl: LUFactorizable where BaseRing: Field & ComputationalRing {
    public func LUfactorize() -> RawLUFactorizationResult {
        let f = _SparseLUFactorizer(self)
        f.debug = true
        f.run()
        return f.result
    }
}

//extension DefaultMatrixImpl: LUFactorizable where BaseRing: ComputationalField {}
