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
    var exclusion: Exclusion
    
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
        self.exclusion = .empty
        
        if exclusion {
            self.exclusion = Exclusion(self)
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
        
        let excDirs = exclusion.excludedDirections
        if v.enumerated().contains(where: { (r, b) in
            excDirs.contains(r) && b == 0
        }) {
            return .zeroModule
        }
        
        let remain = (0 ..< dim).subtract(exclusion.excludedIndeterminates)
        let mons = KR.EdgeRing<R>.monomials(
            ofDegree: 2 * deg,
            usingIndeterminates: remain
        ).map {
            BaseModule.Generator(exponent: $0.leadExponent)
        }
        
        return ModuleStructure<BaseModule>(rawGenerators: mons)
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
            let p = exclusion.exclude(edgeFactor(r))
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
                let f = edge(from: v, to: w)
                print("\t->", w, ":", f(.unit))
            }
        }
    }
    
    //  Simplify chain complex by iterating retractions:
    //             f
    //    Cone( R ---> R )  ==>  Cone( 0 ---> R/f )
    //
    
    struct Exclusion {
        typealias Cube = KRHorizontalCube<R>
        
        private typealias Step = (
            direction: Int,
            excluded: Int,
            divisor: KR.EdgeRing<R>
        )
        
        private var path: [Step]
        private var table: [Int: KR.EdgeRing<R>]
        private var factors: [[Int: KR.EdgeRing<R>]]
        
        private init(_ path: [Step], _ table: [Int: KR.EdgeRing<R>], _ factors: [[Int: KR.EdgeRing<R>]]) {
            self.path = path
            self.table = table
            self.factors = factors
        }
        
        init(_ cube: Cube) {
            typealias P = KR.EdgeRing<R>
            
            var path: [Step] = []
            var table: [Int: KR.EdgeRing<R>] = .empty
            var factors: [[Int: KR.EdgeRing<R>]] = []
            
            let n = cube.dim
            var current = Dictionary(keys: 0 ..< n) { r in
                cube.edgeFactor(r)
            }
            
            factors.append(current)
            
            for r in 0 ..< n {
                let f = current[r]!
                if let i = f.indeterminateOfPrimaryExclusion {

                    table[i] = .indeterminate(i)
                    table = table.mapValues{ $0.divide(by: f, as: i).remainder }
                    
                    current[r] = nil
                    current = current.mapValues { $0.divide(by: f, as: i).remainder }
                    
                    path.append((
                        direction: r,
                        excluded: i,
                        divisor: f
                    ))
                    
                    factors.append(current)
                }
            }
            
            self.init(path, table, factors)
        }
        
        static var empty: Self {
            .init([], [:], [])
        }
        
        var excludedDirections: Set<Int> {
            Set(path.map{ $0.direction })
        }
        
        var excludedIndeterminates: Set<Int> {
            Set(path.map{ $0.excluded })
        }
        
        func exclude(_ f: KR.EdgeRing<R>) -> KR.EdgeRing<R> {
            table.isEmpty ? f : f.substitute(table)
        }
        
        // z ↦ φ(z)
        func collapse(_ z: IndexedModule<Coords, BaseModule>) -> IndexedModule<Coords, BaseModule> {
            if path.isEmpty {
                return z
            }
            
            typealias P = KR.EdgeRing<R>
            
            let excDirs = excludedDirections
            return z.filter{ (v, x) in
                !v.enumerated().contains(where: { (r, b) in
                    excDirs.contains(r) && b == 0
                })
            }.mapValues { x in
                exclude(P(x)).asLinearCombination
            }
        }
        
        // z ↦ ψ(z), chain homotopy inverse of φ.
        func expand(cycle: IndexedModule<Coords, BaseModule>) -> IndexedModule<Coords, BaseModule> {
            if path.isEmpty {
                return cycle
            }
            
            typealias M = IndexedModule<Coords, BaseModule>
            typealias P = KR.EdgeRing<R>
            
            var step = path.count - 1
            var z = Dictionary(keys: cycle.elements.keys) { v in
                P(cycle.elements[v]!)
            }

            for (r, i, f) in path.reversed() {
                // one step back in direction: r.
                let z0 = z.mapPairs { (v, x) -> (Coords, P) in
                    let u = v.replaced(with: 0, at: r)
                    let e = Cube.edgeSign(from: u, to: v)
                    return (u, e * x)
                }

                // w = d'(z0) / f.
                let w = differentiate(z0, step).mapValues { x in
                    x.divide(by: f, as: i).quotient
                }
                
                z = z.merging(w, uniquingKeysWith: +)
                step -= 1
            }
            
            return .init(elements: z.mapValues{ $0.asLinearCombination })
        }
        
        private func differentiate(_ z: [Coords : KR.EdgeRing<R>], _ step: Int) -> [Coords : KR.EdgeRing<R>] {
            typealias P = KR.EdgeRing<R>
            
            let r0 = path[step].direction
            let factors = self.factors[step]
            
            let elements = z.flatMap { (v, x) -> [(Coords, KR.EdgeRing<R>)] in
                v.enumerated().compactMap { (r, b) in
                    if r == r0 || !(factors.contains(key: r) && b == 0) {
                        return nil
                    }
                    let w = v.replaced(with: 1, at: r)
                    let e = Cube.edgeSign(from: v, to: w)
                    let f = factors[r]!
                    let y = e * f * x
                    return (w, y)
                }
            }
            return Dictionary(elements, uniquingKeysWith: +)
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
                generators: Hi.generators.map{ exclusion.expand(cycle: $0) },
                vectorizer: { z in
                    Hi.vectorize( exclusion.collapse(z) )
                }
            )
        }
    }
}
