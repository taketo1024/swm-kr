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
        excluded: Int,
        divisor: KR.EdgeRing<R>
    )]
    private var exclusionTable: [Int: KR.EdgeRing<R>]
    private var edgeFactors: [[Int: KR.EdgeRing<R>]]
    
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
        self.exclusionTable = [:]
        self.edgeFactors = []
        
        if exclusion {
            self.excludeIndeterminates()
        } else {
            // TODO
            edgeFactors = [Dictionary(keys: 0 ..< dim) { r in
                edgeFactor(r)
            }]
        }
    }
    
    //           f
    //  Cone( R ---> R )  ==>  Cone( 0 ---> R/f )
    
    mutating func excludeIndeterminates() {
        typealias P = KR.EdgeRing<R>
        
        let n = L.crossingNumber
        var factors = Dictionary(keys: 0 ..< dim) { r in
            edgeFactor(r)
        }
        var table: [Int: KR.EdgeRing<R>] = .empty
        
        edgeFactors.append(factors)
        
        for r in 0 ..< n {
            let f = factors[r]!
            if f.isLinear && f.degree == 2 { // recall: each xi has deg = 2.
                
                // f = a x_i + (terms) ~ 0
                // <=> x_i ~ -a^-1 (terms)
                
                let a = f.leadCoeff
                let i = f.leadTerm.indexOfIndeterminate
                let y = -a.inverse! * (f - f.leadTerm)
                
//                print("exclude: \(P.indeterminate(i)), f = \(f)")
                
                table = table.mapValues{ $0.substitute([i : y]).reduced }
                table[i] = y
                
                factors[r] = nil
                factors = factors.mapValues { $0.substitute([i : y]) }
                
                exclusions.append((
                    direction: r,
                    excluded: i,
                    divisor: f
                ))
                
                edgeFactors.append(factors)
            }
        }
        
        self.exclusionTable = table
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
        
        let excDirs = excludedDirections
        if v.enumerated().contains(where: { (r, b) in
            excDirs.contains(r) && b == 0
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
        Set(exclusions.map{ $0.excluded })
    }
    
    private func exclude(_ f: KR.EdgeRing<R>) -> KR.EdgeRing<R> {
        f.substitute(exclusionTable)
    }
    
    private func edgeFactor(_ r: Int) -> KR.EdgeRing<R> {
        let c = connection[r]!
        let (ik, il) = (c.ik, c.il)
        
        switch (L.crossings[r].crossingSign, vCoords[r]) {
        case (+1, 0), (-1, 1):
            return ik * il
        case (+1, 1), (-1, 0):
            return ik
        default:
            fatalError("impossible")
        }
    }
    
    func edge(from: Coords, to: Coords) -> ModuleEnd<BaseModule> {
        assert((to - from).weight == 1)
        return edgeCache.getOrSet(key: from.concat(with: to)) {
            let e = edgeSign(from: from, to: to)
            let r = (to - from).enumerated().first { (_, b) in b == 1 }!.offset
            let p = exclude(edgeFactor(r))
            return .init { z -> BaseModule in
                let q = MultivariatePolynomial(z)
                return e * (p * q).asLinearCombination
            }
        }
    }
    
    // z ↦ φ(z)
    private func collapse(cycle: IndexedModule<Coords, BaseModule>) -> IndexedModule<Coords, BaseModule> {
        if exclusions.isEmpty {
            return cycle
        }
        
        typealias P = KR.EdgeRing<R>
        
        let excDirs = excludedDirections
        return cycle.filter{ (v, x) in
            !v.enumerated().contains(where: { (r, b) in
                excDirs.contains(r) && b == 0
            })
        }.mapValues { x in
            exclude(P(x)).asLinearCombination
        }
    }
    
    // z ↦ ψ(z), chain homotopy inverse of φ.
    private func expand(cycle: IndexedModule<Coords, BaseModule>) -> IndexedModule<Coords, BaseModule> {
        if exclusions.isEmpty {
            return cycle
        }
        
        typealias M = IndexedModule<Coords, BaseModule>
        typealias P = KR.EdgeRing<R>
        
        var z = cycle
        var step = exclusions.count - 1
        var dirs = (0 ..< dim).subtract(excludedDirections)

//        print("compute psi of \(z)")
//        print("exclusions: \(exclusions.map{ (d, f, i, _) in (d, i, f) })")
//        print()

        for (r, i, f) in exclusions.reversed() {
            // one step back in direction: r.
            let z0 = z.map { (v, x) in
                let u = v.replaced(with: 0, at: r)
                let e = edgeSign(from: u, to: v)
                return (u, e * x)
            }

            // w = d'(z0) / f.
            let w = differentiate(
                z0,
                step: step,
                movableDirs: dirs
            ).mapValues{
                P($0).divide(by: f, as: i).asLinearCombination
            }
            
            z = (z + w).reduced

//            print("step: \(step), r: \(r), i: \(i), f: \(f)")
//            print("z = \(z)")
//            print("\t-> \(z0)")
//            print("\t-> \(w)")
//            print("z = \(z)\n")

            step -= 1
            dirs.append(r)
        }
        
//        assert(collapse(cycle: z) == cycle)
//        assert(differentiate(z, step: -1, movableDirs: (0 ..< dim).toArray()).isZero)
        
        return z
    }
    
    private func differentiate(_ z: IndexedModule<Coords, BaseModule>, step: Int, movableDirs: [Int]) -> IndexedModule<Coords, BaseModule> {
        typealias P = KR.EdgeRing<R>
        
        let factors = edgeFactors[step]
        return z.elements.sum { (v, x) in
            v.enumerated().filter { (i, b) in
                movableDirs.contains(i) && b == 0
            }.sum { (r, b) -> IndexedModule<Coords, BaseModule> in
                let w = v.replaced(with: 1, at: r)
                let e = edgeSign(from: v, to: w)
                let f = factors[r]!
                let y = e * f * P(x)
                return .init(index: w, value: y.asLinearCombination)
            }
        }
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

extension KRHorizontalCube where R: HomologyCalculatable {
    func homology() -> IndexedModuleStructure<Int, KR.HorizontalModule<R>> {
        let C = self.asChainComplex()
        let H = C.homology()
        
        return .init(support: H.support) { i -> ModuleStructure<KR.HorizontalModule<R>> in
            let Hi = H[i]
            return .init(
                generators: Hi.generators.map{ expand(cycle: $0) },
                vectorizer: { z in
                    Hi.vectorize( collapse(cycle: z) )
                }
            )
        }
    }
}
