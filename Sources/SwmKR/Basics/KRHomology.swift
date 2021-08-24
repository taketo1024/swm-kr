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
    public let exclusion: Bool
    public let symmetry: Bool
    
    private var connection: [Int : KR.EdgeConnection<R>]
    
    private let horizontalHomologyCache: Cache<hKey, HorizontalHomology> = .empty
    private let verticalHomologyCache: Cache<vKey, VerticalHomology> = .empty

    public init(_ L: Link, exclusion: Bool = true, symmetry: Bool = true) {
        self.L = L
        self.exclusion = exclusion
        self.symmetry = symmetry
        self.connection = KREdgeConnection(L).compute()
    }
    
    public subscript(idx: Index) -> Object {
        let (i, j, k) = idx.triple
        if symmetry {
            if !iRange.contains(i) || !jRange.contains(j) || !kRange.contains(k) || !kRange.contains(k + 2 * i) {
                return .zeroModule
            }
            if i > 0 {
                return self[-i, j, k + 2 * i]
            }
        }
        guard let (h, v, s) = ijk2hvs(i, j, k) else {
            return .zeroModule
        }
        let H = verticalHomology(hDegree: h, slice: s)
        return H[v]
    }
    
    public func ijk2hvs(_ i: Int, _ j: Int, _ k: Int) -> (Int, Int, Int)? {
        let g = baseGrading
        let (a, b, c) = (g[0], g[1], g[2])

        if (i - a).isEven && (j - b).isEven && (k - c).isEven {
            return ((j - b)/2, (k - c)/2, (i - a)/2 - (j - b)/2)
        } else {
            return nil
        }
    }
    
    public func hvs2ijk(_ h: Int, _ v: Int, _ s: Int) -> (Int, Int, Int) {
        let g = baseGrading
        let (a, b, c) = (g[0], g[1], g[2])
        return (a + 2 * h + 2 * s,
                b + 2 * h,
                c + 2 * v)
    }
    
    public var baseGrading: KR.Grading {
        let n = L.crossingNumber
        let w = L.writhe
        let s = L.numberOfSeifertCircles

        let v0 = Cube.Coords.zeros(length: n)
        let base = KR.baseGrading(link: L, hCoords: v0, vCoords: v0)
        let gradingShift: KR.Grading = [-w + s - 1, w + s - 1, w - s + 1]
        
        return base + gradingShift
    }
    
    public var levelRange: ClosedRange<Int> {
        let n = L.crossingNumber
        let s = L.numberOfSeifertCircles
        
        let l0 = -2 * n
        let l1 = (-3 * n + s - 1) / 2

        return l0 ... l1
    }
    
    public var iRange: ClosedRange<Int> {
        let n = L.crossingNumber
        let s = L.numberOfSeifertCircles
        
        let i0 = -n + s - 1
        let i1 =  n - s + 1
        
        return i0 ... i1
    }
    
    public var jRange: ClosedRange<Int> {
        let w = L.writhe
        let s = L.numberOfSeifertCircles

        let j0 = w - s + 1
        let j1 = w + s - 1
        
        return j0 ... j1
    }
    
    public var kRange: ClosedRange<Int> {
        let n = L.crossingNumber
        
        let k0 = baseGrading[2]
        let k1 = k0 + 2 * n
        
        return k0 ... k1
    }

    public var support: [MultiIndex<_3>] {
        let n = L.crossingNumber
        return levelRange.flatMap { s in
            ((0 ... n) * (0 ... n)).map { (h, v) in
                let (i, j, k) = hvs2ijk(h, v, s)
                return [i, j, k]
            }
        }
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
            return C.homology(options: .onlyStructures)
        }
    }
    
    public func clearCache() {
        horizontalHomologyCache.clear()
        verticalHomologyCache.clear()
    }
}

// convenience methods
extension KRHomology {
    public func structure(restrictedTo: (Int, Int, Int) -> Bool = { (_, _, _) in true }) -> [KR.Grading : ModuleStructure<KR.TotalModule<R>>] {
        typealias E = (KR.Grading, ModuleStructure<KR.TotalModule<R>>)
        
        let n = L.crossingNumber
        let elements = levelRange.flatMap { s -> [E] in
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
    
    public var lowestQPart: KR.qaPolynomial<𝐙> {
        let str = structure { (i, _, _) in i == iRange.lowerBound }
        return qaPolynomial(str)
    }
    
    public var highestQPart: KR.qaPolynomial<𝐙> {
        let str = structure { (i, _, _) in i == iRange.upperBound }
        return qaPolynomial(str)
    }
    
    public var lowestAPart: KR.qaPolynomial<𝐙> {
        let str = structure { (_, j, _) in j == jRange.lowerBound }
        return qaPolynomial(str)
    }
    
    public var highestAPart: KR.qaPolynomial<𝐙> {
        let str = structure { (_, j, _) in j == jRange.upperBound }
        return qaPolynomial(str)
    }
    
    private func qaPolynomial(_ structure: [KR.Grading : ModuleStructure<KR.TotalModule<R>>]) -> KR.qaPolynomial<𝐙> {
        .init(elements: structure.map { (g, V) in
            let (i, j, k) = (g[0], g[1], g[2])
            return ([i, j], (-1).pow( (k - j) / 2) * V.rank )
        })
    }
}

fileprivate struct hKey: Hashable {
    let vCoords: Cube.Coords
    let slice: Int
}

fileprivate struct vKey: Hashable {
    let hDegree: Int
    let slice: Int
}
