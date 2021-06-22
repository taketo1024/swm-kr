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
    
    private var exclusions: [(
        direction: Int,
        factor: KR.EdgeRing<R>,
        exclude: Int,
        table: [Int: KR.EdgeRing<R>]
    )]
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
        
        self.exclusions = []
        
        if exclusion {
            self.excludeIndeterminates()
        }
    }
    
    //           f
    //  Cone( R ---> R )  ==>  Cone( 0 ---> R/f )
    
    mutating func excludeIndeterminates() {
        typealias P = KR.EdgeRing<R>
        
        let n = L.crossingNumber
        var table: [Int: KR.EdgeRing<R>] = .empty
        
        for c in 0 ..< n {
            let f = edgeFactor(c)
            if f.isLinear && f.degree == 2 { // recall: each xi has deg = 2.
                
                // f = a x_i + (terms) ~ 0
                // <=> x_i ~ -a^-1 (terms)
                
                let a = f.leadCoeff
                let i = f.leadTerm.indexOfIndeterminate
                let y = -a.inverse! * (f - f.leadTerm)
                
//                print("dir: \(c), exclude: \(P.indeterminate(i)), f: \(f)")
                
                table = table.mapValues{ $0.substitute([i : y]).reduced }
                table[i] = y
                
                exclusions.append((
                    direction: c,
                    factor: f,
                    exclude: i,
                    table: table
                ))
            }
        }
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
        
        if v.enumerated().contains(where: { (i, b) in
            exclusions.contains(where: { $0.direction == i }) && b == 0
        }) {
            return .zeroModule
        }
        
        
        let remain = (0 ..< dim).subtract(excludedIndeterminates)
        let mons = KR.EdgeRing<R>.monomials(
            ofDegree: 2 * deg,
            usingIndeterminates: remain
        ).map {
            BaseModule.Generator(exponent: $0.leadExponent)
        }
        
        return ModuleStructure<BaseModule>(rawGenerators: mons)
    }
    
    private var excludedDirections: Set<Int> {
        Set(exclusions.map{ $0.direction })
    }
    
    private var excludedIndeterminates: Set<Int> {
        if let indices = exclusions.last?.table.keys {
            return Set(indices)
        } else {
            return []
        }
    }
    
    private func exclude(_ f: KR.EdgeRing<R>, step: Int? = nil) -> KR.EdgeRing<R> {
        let table = step.map{ i in i >= 0 ? exclusions[i].table : [:] }
            ?? exclusions.last?.table
            ?? [:]
        
        return f.substitute(table)
    }
    
    private func edgeFactor(_ i: Int, step: Int? = nil) -> KR.EdgeRing<R> {
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

        return exclude(f, step: step)
    }
    
    func edge(from: Coords, to: Coords) -> ModuleEnd<BaseModule> {
        assert((to - from).weight == 1)
        return edgeCache.getOrSet(key: from.concat(with: to)) {
            let e = edgeSign(from: from, to: to)
            let i = (to - from).enumerated().first { (_, b) in b == 1 }!.offset
            let p = edgeFactor(i)
            return .init { z -> BaseModule in
                let q = MultivariatePolynomial(z)
                return e * (p * q).asLinearCombination
            }
        }
    }
    
    func recover(cycle: IndexedModule<Coords, BaseModule>) -> IndexedModule<Coords, BaseModule> {
        // TODO
        .zero
    }
    
    func printDescription() {
        for v in BitSequence.allSequences(length: dim) {
            print(v, ":", self[v].generators)
            for w in v.successors {
                let f = edge(from: v, to: w)
                print("\t->", w, ":", f(.unit))
            }
        }
    }
}
