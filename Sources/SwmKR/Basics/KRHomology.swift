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

public struct KRHomology<R: HomologyCalculatable>: GradedModuleStructureType {
    public typealias BaseModule = KR.TotalModule<R>
    public typealias Index = IntList<_3>
    public typealias Object = ModuleStructure<BaseModule>

    public typealias HorizontalHomology = GradedModuleStructure<Int, KR.HorizontalModule<R>>
    public typealias VerticalHomology = GradedModuleStructure<Int, KR.TotalModule<R>>

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
        let gradingShift = KR.Grading(-w + s - 1, w + s - 1, w - s + 1)
        
        return base + gradingShift
    }
    
    public var levelRange: ClosedRange<Int> {
        let n = L.crossingNumber
        return -2 * n ... -n
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
        iRange
    }

    public var support: [IntList<_3>] {
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

fileprivate struct hKey: Hashable {
    let vCoords: Cube.Coords
    let slice: Int
}

fileprivate struct vKey: Hashable {
    let hDegree: Int
    let slice: Int
}
