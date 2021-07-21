//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/05/10.
//

import SwmCore
import SwmKnots
import SwmHomology
import SwmKhovanov

public struct KRHomology<R: HomologyCalculatable>: IndexedModuleStructureType {
    public typealias BaseModule = KR.TotalModule<R>
    public typealias Index = MultiIndex<_3>
    public typealias Object = ModuleStructure<BaseModule>

    public typealias HorizontalHomology = IndexedModuleStructure<Int, KR.HorizontalModule<R>>
    public typealias VerticalHomology = IndexedModuleStructure<Int, KR.TotalModule<R>>

    public let L: Link
    public let normalized: Bool
    public let exclusion: Bool
    public let symmetry: Bool
    
    private let gradingShift: KR.Grading
    private var connection: [Int : KR.EdgeConnection<R>]
    
    private let horizontalHomologyCache: Cache<hKey, HorizontalHomology> = .empty
    private let verticalHomologyCache: Cache<vKey, VerticalHomology> = .empty

    public init(_ L: Link, normalized: Bool = true, exclusion: Bool = true, symmetry: Bool = true) {
        let w = L.writhe
        let b = L.resolved(by: L.orientationPreservingState).components.count
        
        self.L = L
        self.normalized = normalized
        self.exclusion = exclusion
        self.symmetry = symmetry
        self.gradingShift = normalized ? [-w + b - 1, w + b - 1, w - b + 1] : .zero
        self.connection = KREdgeConnection(L).compute()
    }
    
    public subscript(idx: Index) -> Object {
        let (i, j, k) = idx.triple
        if normalized && symmetry {
            if i < lowestQ || highestQ < i {
                return .zeroModule
            }
            if j < lowestA || highestA < j {
                return .zeroModule
            }
            if i < 0 {
                return self[-i, j, k + 2 * i]
            }
        }
        guard let (h, v, s) = ijk2hvs(i, j, k) else {
            return .zeroModule
        }
        let H = verticalHomology(hDegree: h, slice: s)
        return H[v]
    }
    
    private var minSlice: Int {
        -2 * L.crossingNumber
    }
    
    private var highestQ: Int {
        L.crossingNumber - L.numberOfSeifertCircles + 1
    }
    
    private var lowestQ: Int {
        -L.crossingNumber + L.numberOfSeifertCircles - 1
    }
    
    private var highestA: Int {
        L.writhe + L.numberOfSeifertCircles - 1
    }
    
    private var lowestA: Int {
        L.writhe - L.numberOfSeifertCircles + 1
    }
    
    public var support: [MultiIndex<_3>] {
        let n = L.crossingNumber
        return ((0 ... n) * (0 ... n)).flatMap { (h, v) in
            (minSlice ... 0).map { s in
                let (i, j, k) = hvs2ijk(h, v, s)
                return [i, j, k]
            }
        }
    }
    
    public var baseGrading: KR.Grading {
        let n = L.crossingNumber
        let v0 = Cube.Coords.zeros(length: n)
        return KR.baseGrading(link: L, hCoords: v0, vCoords: v0) + gradingShift
    }
    
    public func grading(of z: KR.TotalModule<R>) -> KR.Grading {
        if z.isZero {
            return .zero
        }
        
        let (v, x) = z.elements.anyElement!
        let (h, y) = x.elements.anyElement!
        let q = y.degree
        
        let g = KR.baseGrading(link: L, hCoords: h, vCoords: v)
        return g + [2 * q, 0, 0] + gradingShift
    }
    
    public func horizontalComplex(at vCoords: Cube.Coords, slice: Int) -> ChainComplex1<KR.HorizontalModule<R>> {
        let cube = KRHorizontalCube(link: L, vCoords: vCoords, slice: slice, connection: connection, exclusion: exclusion)
        return cube.asChainComplex()
    }
    
    public func horizontalHomology(at vCoords: Cube.Coords, slice: Int) -> HorizontalHomology {
        let key = hKey(vCoords: vCoords, slice: slice)
        return horizontalHomologyCache.getOrSet(key: key) {
            let cube = KRHorizontalCube(link: L, vCoords: vCoords, slice: slice, connection: connection, exclusion: exclusion)
            return cube.homology()
        }
    }
    
    public func verticalComplex(hDegree h: Int, slice s: Int) -> ChainComplex1<KR.TotalModule<R>> {
        let cube = KRTotalCube<R>(link: L, connection: connection) { vCoords -> KRTotalCube<R>.Vertex in
            let H = horizontalHomology(at: vCoords, slice: s)
            if !H[h].isFree {
                fatalError("Horizontal homology at \(vCoords) is non-free.")
            }
            return H[h]
        }
        return cube.asChainComplex()
    }
    
    public func verticalHomology(hDegree h: Int, slice s: Int) -> VerticalHomology {
        let key = vKey(hDegree: h, slice: s)
        return verticalHomologyCache.getOrSet(key: key) {
            let C = verticalComplex(hDegree: h, slice: s)
            return C.homology()
        }
    }
    
    public func structure(restrictedTo: (Int, Int, Int) -> Bool = { (_, _, _) in true }) -> [KR.Grading : ModuleStructure<KR.TotalModule<R>>] {
        typealias E = (KR.Grading, ModuleStructure<KR.TotalModule<R>>)
        
        let n = L.crossingNumber
        let elements = (minSlice ... 0).flatMap { s -> [E] in
            ((0 ... n) * (0 ... n)).compactMap { (h, v) -> E? in
                let (i, j, k) = hvs2ijk(h, v, s)
                if restrictedTo(i, j, k) {
                    let H = self[i, j, k]
                    return !H.isZero ? ([i, j, k], H) : nil
                } else {
                    return nil
                }
            }
        }
        return Dictionary(elements)
    }
    
    public func table(restrictedTo: (Int, Int, Int) -> Bool = { (_, _, _) in true }) -> String {
        typealias qPoly = KR.qPolynomial<𝐙>
        let str = structure(restrictedTo: restrictedTo)
        let table = str
            .group { (ijk, _) in
                MultiIndex<_2>(ijk[1], ijk[2])
            }
            .mapValues{ list in
                qPoly(elements: list.map{ (ijk, V) in
                    (ijk[0], V.rank)
                })
            }
        
        return Format.table(
            rows: table.keys.map{ $0[1] }.uniqued().sorted().reversed(), // k
            cols: table.keys.map{ $0[0] }.uniqued().sorted(), // j
            symbol: "k\\j",
            printHeaders: true
        ) { (k, j) -> String in
            let q = table[[j, k]] ?? .zero
            return !q.isZero ? "\(q)" : ""
        }
    }
    
    public func table(i: Int? = nil, j: Int? = nil, k: Int? = nil) -> String {
        table { (i1, j1, k1) in
            return
                (i.flatMap{ $0 == i1 } ?? true)
             && (j.flatMap{ $0 == j1 } ?? true)
             && (k.flatMap{ $0 == k1 } ?? true)
        }
    }
    
    public func printTable(restrictedTo r: (Int, Int, Int) -> Bool = { (_, _, _) in true }) {
        print(table(restrictedTo: r))
    }
    
    public func printTable(i: Int? = nil, j: Int? = nil, k: Int? = nil) {
        print(table(i: i, j: j, k: k))
    }
    
    public var gradedEulerCharacteristic: KR.qaPolynomial<𝐙> {
        qaPolynomial(structure())
    }
    
    public var highestQPart: KR.qaPolynomial<𝐙> {
        let str = structure { (i, _, _) in i == highestQ }
        return qaPolynomial(str)
    }
    
    public var lowestQPart: KR.qaPolynomial<𝐙> {
        let str = structure { (i, _, _) in i == lowestQ }
        return qaPolynomial(str)
    }
    
    public var highestAPart: KR.qaPolynomial<𝐙> {
        let str = structure { (_, j, _) in j == highestA }
        return qaPolynomial(str)
    }
    
    public var lowestAPart: KR.qaPolynomial<𝐙> {
        let str = structure { (_, j, _) in j == lowestA }
        return qaPolynomial(str)
    }
    
    private func qaPolynomial(_ structure: [KR.Grading : ModuleStructure<KR.TotalModule<R>>]) -> KR.qaPolynomial<𝐙> {
        .init(elements: structure.map { (g, V) in
            let (i, j, k) = (g[0], g[1], g[2])
            return ([i, j], (-1).pow( (k - j) / 2) * V.rank )
        })
    }
    
    private func ijk2hvs(_ i: Int, _ j: Int, _ k: Int) -> (Int, Int, Int)? {
        let g = baseGrading
        let (a, b, c) = (g[0], g[1], g[2])

        if (i - a).isEven && (j - b).isEven && (k - c).isEven {
            return ((j - b)/2, (k - c)/2, (i - a)/2 - (j - b)/2)
        } else {
            return nil
        }
    }
    
    private func hvs2ijk(_ h: Int, _ v: Int, _ s: Int) -> (Int, Int, Int) {
        let g = baseGrading
        let (a, b, c) = (g[0], g[1], g[2])
        return (a + 2 * h + 2 * s,
                b + 2 * h,
                c + 2 * v)
    }
    
    private struct hKey: Hashable {
        let vCoords: Cube.Coords
        let slice: Int
    }

    private struct vKey: Hashable {
        let hDegree: Int
        let slice: Int
    }
}
