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

internal struct KRTotalCube<R: Ring>: ModuleCube {
    typealias BaseModule = KR.HorizontalModule<R>
    typealias Vertex = ModuleStructure<BaseModule>
    typealias Edge = ModuleEnd<BaseModule>

    let L: Link
    let connection: [Int : KR.EdgeConnection<R>]
    let vertex: (Coords) -> Vertex
    
    init(link L: Link, connection: [Int : KR.EdgeConnection<R>], vertex: @escaping (Coords) -> Vertex) {
        self.L = L
        self.connection = connection
        self.vertex = vertex
    }
    
    var dim: Int {
        L.crossingNumber
    }
    
    subscript(v: Coords) -> ModuleStructure<BaseModule> {
        vertex(v)
    }
    
    private func edgeFactor(_ p: Int, subcoords: Coords) -> KR.EdgeRing<R> {
        let il = connection[p]!.il
        
        switch (L.crossings[p].crossingSign, subcoords[p]) {
        case (+1, 0), (-1, 1):
            return il
        case (+1, 1), (-1, 0):
            return .identity
        default:
            fatalError("impossible")
        }
    }
    
    func edge(from: Coords, to: Coords) -> ModuleEnd<BaseModule> {
        let e = edgeSign(from: from, to: to)
        let i = (to - from).enumerated().first { (_, b) in b == 1 }!.offset
        return .init { x -> BaseModule in
            x.elements.sum { (subcoords, z) -> BaseModule in
                let p = edgeFactor(i, subcoords: subcoords)
                let q = MultivariatePolynomial(z)
                return IndexedModule(
                    index: subcoords,
                    value: e * (p * q).asLinearCombination
                )
            }
        }
    }
}
