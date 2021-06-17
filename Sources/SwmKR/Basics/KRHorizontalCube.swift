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

internal struct KRHorizontalCube<R: Ring>: ModuleCube {
    typealias BaseModule = KR.BaseModule<R>
    typealias Vertex = ModuleStructure<BaseModule>
    typealias Edge = ModuleEnd<BaseModule>

    let L: Link
    let vCoords: Coords // Cubic coord in Cube2
    let slice: Int
    let connection: [Int : KR.EdgeConnection<R>]
    
    private var excludedDirections: Set<Int>
    private var excludedIndeterminates: [Int: KR.EdgeRing<R>]
    
    private let vertexCache: Cache<Coords, Vertex> = .empty
    private let   edgeCache: Cache<Coords, Edge>   = .empty
    
    init(link L: Link, vCoords: Coords, slice: Int, exclusion: Bool = false) { // TODO 
        let conn = KREdgeConnection<R>(L).compute()
        self.init(link: L, vCoords: vCoords, slice: slice, connection: conn, exclusion: exclusion)
    }
    
    init(link L: Link, vCoords: Coords, slice: Int, connection: [Int : KR.EdgeConnection<R>], exclusion: Bool = false) {
        self.L = L
        self.vCoords = vCoords
        self.slice = slice
        self.connection = connection
        
        self.excludedDirections = []
        self.excludedIndeterminates = [:]
        
        if exclusion {
            self.exclusion()
        }
    }
    
    //           f
    //  Cone( R ---> R )  ==>  Cone( 0 ---> R/f )
    
    mutating func exclusion() {
        typealias P = KR.EdgeRing<R>
        
        let n = L.crossingNumber
        for c in 0 ..< n {
            let f = edgeFactor(c)
            if f.isLinear && f.degree == 2 {
                excludedDirections.insert(c)
                
                // f = a x_i + (terms) ~ 0
                // <=> x_i ~ -a^-1 (terms)
                
                let a = f.leadCoeff
                let i = f.leadTerm.indexOfIndeterminate
                let y = -a.inverse! * (f - f.leadTerm)
                
                excludedIndeterminates = excludedIndeterminates.mapValues {
                    $0.substitute([i: y])
                }
                excludedIndeterminates[i] = y
                
//                print("exclude \(P.indeterminate(i)) -> \(y)")
            }
        }
//        print("excludeded-directions:", excludedDirections)
    }
    
    var dim: Int {
        L.crossingNumber
    }
    
    var baseGrading: KR.Grading {
        let v0 = Coords.zeros(length: dim)
        return KR.baseGrading(link: L, hCoords: v0, vCoords: v0)
    }
    
    func gradingShift(at hCoords: Coords) -> KR.Grading {
        KR.baseGrading(link: L, hCoords: hCoords, vCoords: vCoords)
    }
    
    subscript(v: Coords) -> Vertex {
        vertexCache.getOrSet(key: v) {
            vertex(v)
        }
    }
    
    private func vertex(_ v: Coords) -> Vertex {
        let deg = slice + v.weight + (baseGrading - gradingShift(at: v))[0] / 2
        if deg < 0 {
            return .zeroModule
        }
        
        let excluded = v.enumerated().contains{ (i, b) in
            excludedDirections.contains(i) && b == 0
        }
        if excluded {
            return .zeroModule
        }
        
        let indeterminates = (0 ..< dim).subtract(Set(excludedIndeterminates.keys))
        let mons = KR.EdgeRing<R>.monomials(
            ofDegree: 2 * deg,
            usingIndeterminates: indeterminates
        ).map {
            BaseModule.Generator(exponent: $0.leadExponent)
        }
        
        return ModuleStructure<BaseModule>(rawGenerators: mons)
    }
    
    private func edgeFactor(_ i: Int) -> KR.EdgeRing<R> {
        let c = connection[i]!
        let (ik, il) = (c.ik, c.il)
        
        let f = { () -> KR.EdgeRing<R> in
            switch (L.crossings[i].crossingSign, vCoords[i]) {
            case (+1, 0), (-1, 1):
                return ik * il
            case (+1, 1), (-1, 0):
                return ik
            default:
                fatalError("impossible")
            }
        }()

        return f.substitute(excludedIndeterminates)
    }
    
    private func edgeFactor(from: Coords, to: Coords) -> KR.EdgeRing<R> {
        assert((to - from).weight == 1)
        let i = (to - from).enumerated().first { (_, b) in b == 1 }!.offset
        return edgeFactor(i)
    }
    
    func edge(from: Coords, to: Coords) -> ModuleEnd<BaseModule> {
        edgeCache.getOrSet(key: from.concat(with: to)) {
            let e = edgeSign(from: from, to: to)
            let p = edgeFactor(from: from, to: to)
            return .init { z -> BaseModule in
                let q = MultivariatePolynomial(z)
                return e * (p * q).asLinearCombination
            }
        }
    }
    
    func printDescription() {
        for v in BitSequence.allSequences(length: dim) {
            print(v, ":", self[v].generators)
            for w in v.successors {
                print("\t->", w, ":", edge(from: v, to: w).callAsFunction(.unit))
            }
        }
    }
}
