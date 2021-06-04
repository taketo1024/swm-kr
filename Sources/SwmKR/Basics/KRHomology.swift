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

public struct KRHomology<R: EuclideanRing> {
    public let L: Link
    public let normalized: Bool
    
    private let gradingShift: KR.Grading
    private var connection: [Int : KR.EdgeConnection<R>]
    
    private let horizontalHomologyCache: Cache<hKey, ModuleGrid1<KR.HorizontalModule<R>>> = .empty
    private let   verticalHomologyCache: Cache<vKey, ModuleGrid1<KR.TotalModule<R>>> = .empty

    public init(_ L: Link, normalized: Bool = true) {
        let w = L.writhe
        let b = L.resolved(by: L.orientationPreservingState).components.count
        
        self.L = L
        self.normalized = normalized
        self.gradingShift = normalized ? [-w + b - 1, w + b - 1, w - b + 1] : .zero
        self.connection = KREdgeConnection(L).compute()
    }
    
    public subscript(i: Int, j: Int, k: Int) -> ModuleStructure<KR.TotalModule<R>> {
        guard let (h, v, s) = ijk2hvs(i, j, k) else {
            return .zeroModule
        }
        let H = verticalHomology(hDegree: h, slice: s)
        return H[v]
    }
    
    public var minSlice: Int {
        -2 * L.crossingNumber
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
        let cube = KRHorizontalCube(link: L, vCoords: vCoords, slice: slice, connection: connection)
        return cube.asChainComplex()
    }
    
    public func horizontalHomology(at vCoords: Cube.Coords, slice: Int) -> ModuleGrid1<KR.HorizontalModule<R>> {
        let key = hKey(vCoords: vCoords, slice: slice)
        return horizontalHomologyCache.getOrSet(key: key) {
            let C = self.horizontalComplex(at: vCoords, slice: slice)
            return C.homology()
        }
    }
    
    public func verticalComplex(hDegree h: Int, slice s: Int) -> ChainComplex1<KR.TotalModule<R>> {
        let cube = KRTotalCube<R>(link: L, connection: connection) { vCoords -> KRTotalCube<R>.Vertex in
            let H = horizontalHomology(at: vCoords, slice: s)
            return H[h]
        }
        return cube.asChainComplex()
    }
    
    public func verticalHomology(hDegree h: Int, slice s: Int) -> ModuleGrid1<KR.TotalModule<R>> {
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
