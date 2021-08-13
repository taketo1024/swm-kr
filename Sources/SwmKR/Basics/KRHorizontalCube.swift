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
    
    init(link L: Link, vCoords: Coords, slice: Int, exclusion: Bool = false) { // TODO 
        let conn = KREdgeConnection<R>(L).compute()
        self.init(link: L, vCoords: vCoords, slice: slice, connection: conn, exclusion: exclusion)
    }
    
    init(link L: Link, vCoords: Coords, slice: Int, connection: [Int : KR.EdgeConnection<R>], exclusion: Bool = false) {
        self.L = L
        self.vCoords = vCoords
        self.slice = slice
        self.connection = connection
        self.exclusion = .empty // must initialize first
        self.exclusion = Exclusion(self, exclusion: exclusion)
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
        
        let primExc = exclusion.primarilyExcludedIndeterminates
        let secExc = exclusion.secondarilyExcludedIndeterminates
        let remain = (0 ..< dim).subtract(primExc)
        
        let mons = KR.EdgeRing<R>.monomials(
            ofDegree: 2 * deg,
            usingIndeterminates: remain
        ).exclude { m in
            let e = m.leadExponent
            return secExc.contains { i in
                e[i] >= 2
            }
        }
        
        return ModuleStructure<BaseModule>(rawGenerators: mons.map {
            BaseModule.Generator(exponent: $0.leadExponent)
        })
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
        let r = (to - from).enumerated().first { (_, b) in b == 1 }!.offset
        let p = exclusion.edgeFactor(r)
        if p.isZero {
            return .zero
        } else {
            let e = edgeSign(from: from, to: to)
            return .init { z -> BaseModule in
                let q = MultivariatePolynomial(z)
                let r = e * exclusion.exclude(p * q)
                return r.asLinearCombination
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
        private var factors: [[Int: KR.EdgeRing<R>]]
        
        private init(_ path: [Step], _ factors: [[Int: KR.EdgeRing<R>]]) {
            self.path = path
            self.factors = factors
        }
        
        init(_ cube: Cube, exclusion: Bool = false) {
            typealias P = KR.EdgeRing<R>
            
            var path: [Step] = []
            var factors: [[Int: KR.EdgeRing<R>]] = []
            
            let n = cube.dim
            var current = Dictionary(keys: 0 ..< n) { r in
                cube.edgeFactor(r)
            }
            
            var vars = Set(0 ..< n)
            
            factors.append(current)
            
            func append(_ r: Int, _ i: Int, _ f: KR.EdgeRing<R>) {
                current[r] = nil
                current = current.mapValues {
                    $0.divide(by: f, as: i).remainder
                }.exclude{
                    $0.value.isZero
                }
                
                path.append((
                    direction: r,
                    excluded: i,
                    divisor: f
                ))
                
                factors.append(current)
            }
            
            if exclusion {
                for r in 0 ..< n where current[r] != nil {
                    let f = current[r]!
                    if let i = f.primaryExclusive(in: vars) {
                        append(r, i, f)
                        vars.remove(i)
                    }
                }
                
                for r in 0 ..< n where current[r] != nil {
                    let f = current[r]!
                    if let i = f.secondaryExclusive(in: vars) {
//                        print("exclude x\(Format.sub(i)) by \(f)")
                        append(r, i, f)
                        vars.subtract(f.involvedIndeterminates)
                        if vars.isEmpty {
                            break
                        }
                    }
                }
            }
            
            self.init(path, factors)
        }
        
        static var empty: Self {
            .init([], [])
        }
        
        var excludedDirections: Set<Int> {
            Set(path.map{ $0.direction })
        }
        
        var primarilyExcludedIndeterminates: Set<Int> {
            Set(path.compactMap{ (_, i, f) in
                f.degree(as: i) == 1 ? i : nil
            })
        }
        
        var secondarilyExcludedIndeterminates: Set<Int> {
            Set(path.compactMap{ (_, i, f) in
                f.degree(as: i) == 2 ? i : nil
            })
        }
        
        func exclude(_ g: KR.EdgeRing<R>) -> KR.EdgeRing<R> {
            path.reduce(g) { (g, next) in
                let (_, i, f) = next
                return g.divide(by: f, as: i).remainder
            }
        }
        
        func edgeFactor(_ r: Int) -> KR.EdgeRing<R> {
            factors.last?[r] ?? .zero
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
