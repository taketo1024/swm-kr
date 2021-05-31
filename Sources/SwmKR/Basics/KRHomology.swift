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
    private let hHomologyCache: Cache<hKey, ModuleGrid1<KR.HorizontalModule<R>>> = .empty

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
        let C = totalComplex(hDegree: h, slice: s)
        let H = C.homology()
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
    
    public func totalComplex(hDegree h: Int, slice: Int) -> ChainComplex1<KR.TotalModule<R>> {
        let cube = KRTotalCube<R>(link: L, connection: connection) { vCoords -> KRTotalCube<R>.Vertex in
            let H = hHomologyCache.getOrSet(key: hKey(vCoords: vCoords, slice: slice)) {
                let C = self.horizontalComplex(at: vCoords, slice: slice)
                let H = C.homology()
                return H
            }
            return H[h]
        }
        return cube.asChainComplex()
    }
    
    public func structure() -> [KR.Grading : ModuleStructure<KR.TotalModule<R>>] {
        typealias E = (KR.Grading, ModuleStructure<KR.TotalModule<R>>)
        
        let n = L.crossingNumber
        let seq = Array(minSlice ... 0).flatMap { s -> [E] in
            Array(0 ... n).flatMap { h -> [E] in
                let Cv = totalComplex(hDegree: h, slice: s)
                let H = Cv.homology()
                return (0 ... n).compactMap { v -> E? in
                    let (i, j, k) = hvs2ijk(h, v, s)
                    let Hv = H[v]
                    return !Hv.isZero ? ([i, j, k], Hv) : nil
                }
            }
        }
        return Dictionary(seq)
    }
    
    public func gradedEulerCharacteristic() -> String {
        structure().map { (g, obj) -> String in
            let (i, j, k) = (g[0], g[1], g[2])
            return Format.linearCombination( [ ("a\(Format.sup(j))q\(Format.sup(i))", (-1).pow( (k - j) / 2)) ] )
        }.joined(separator: " + ")
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
        return (a + 2 * h + 2 * s, b + 2 * h, c + 2 * v)
    }
    
    private struct hKey: Hashable {
        let vCoords: Cube.Coords
        let slice: Int
    }
}
