//
//  GridComplex.swift
//  SwiftyHomology
//
//  Created by Taketo Sano on 2019/07/04.
//

import SwiftyMath
import SwiftyHomology

public struct _Un: MPolynomialIndeterminate {
    public static let numberOfIndeterminates = Int.max
    public static func degree(_ i: Int) -> Int {
        return -2
    }
    public static func symbol(_ i: Int) -> String {
        return "U\(Format.sub(i))"
    }
}

public typealias GridComplex = ChainComplex1<FreeModule<TensorGenerator<MonomialGenerator<_Un>, GridDiagram.Generator>, 𝐙₂>>

extension GridComplex {
    // GC-tilde. [Book] p.72, Def 4.4.1
    // MEMO: not a knot invariant.
    public static func tilde(_ G: GridDiagram) -> ChainComplex1<FreeModule<GridDiagram.Generator, 𝐙₂>> {
        typealias R = 𝐙₂
        
        let (Os, Xs) = (G.Os, G.Xs)
        return _gridComplex(G) { rect in
            (rect.intersects(Xs) || rect.intersects(Os))
                ? nil : .identity
        }
    }
    
    // GC-hat. [Book] p.80, Def 4.6.12.
    public static func hat(_ G: GridDiagram) -> GridComplex {
        typealias P = MPolynomial<_Un, 𝐙₂>
        
        let (Os, Xs) = (G.Os, G.Xs)
        let O_last = Os.last!
        
        return _gridComplex(G) { rect in
            (rect.intersects(Xs) || rect.contains(O_last))
                ? nil : P.monomial(ofMultiDegree: Os.map { O in rect.contains(O) ? 1 : 0 })
        }.splitMonomials(numberOfIndeterminants: G.gridNumber - 1)
    }
    
    // GC^-. [Book] p.75, Def 4.6.1
    public static func minus(_ G: GridDiagram) -> GridComplex {
        typealias P = MPolynomial<_Un, 𝐙₂>
        
        let (Os, Xs) = (G.Os, G.Xs)
        return _gridComplex(G) { rect in
            rect.intersects(Xs)
                ? nil : P.monomial(ofMultiDegree: Os.map { O in rect.contains(O) ? 1 : 0 })
        }.splitMonomials(numberOfIndeterminants: G.gridNumber)
    }
    
    // MEMO: 𝓖𝓒^-. [Book] p.252, Def 13.2.1
    public static func filtered(_ G: GridDiagram) -> GridComplex {
        typealias P = MPolynomial<_Un, 𝐙₂>
        
        let Os = G.Os
        return _gridComplex(G) { rect in
            P.monomial(ofMultiDegree: Os.map { O in rect.contains(O) ? 1 : 0 })
        }.splitMonomials(numberOfIndeterminants: G.gridNumber)
    }
    
    private static func _gridComplex<R: Ring>(_ G: GridDiagram, _ coeff: @escaping (GridDiagram.Rect) -> R?) -> ChainComplex1<FreeModule<GridDiagram.Generator, R>> {
        return ChainComplex1.descending(
            supported: G.generators.map{ $0.MaslovDegree }.range!,
            sequence: { i in ModuleObject(basis: G.generators.filter{ $0.degree == i }) },
            differential: { i in
                ModuleEnd.linearlyExtend {  x -> FreeModule<GridDiagram.Generator, R> in
                    G.adjacents(x).sum { y -> FreeModule<GridDiagram.Generator, R> in
                        let rects = G.emptyRectangles(from: x, to: y)
                        let c = rects.compactMap { rect in coeff(rect) }.sumAll()
                        return (c == .zero) ? .zero : c * .wrap(y)
                    }
                }
            }
        )
    }
}

extension ChainComplex where GridDim == _1, BaseModule == FreeModule<GridDiagram.Generator, MPolynomial<_Un, 𝐙₂>> {
    func splitMonomials(numberOfIndeterminants n: Int) -> GridComplex {
        typealias R = 𝐙₂
        typealias P = MPolynomial<_Un, R>
        typealias Result = FreeModule<TensorGenerator<MonomialGenerator<_Un>, GridDiagram.Generator>, R>
        
        let iMax = grid.supportedCoords.map{ $0[0] }.max()!
        return ChainComplex1<Result>.descending(
            supported: grid.supportedCoords.map{ $0[0] },
            sequence: { i -> ModuleObject<Result> in
                guard i <= iMax else {
                    return .zeroModule
                }
                
                let above = (0 ... (iMax - i) / 2).flatMap { k in self[i + 2 * k].generators }
                let gens = above.flatMap { e -> [Result.Generator] in
                    let x = e.unwrap()
                    let mons = P.monomials(ofDegree: i - x.degree, usingIndeterminates: (0 ..< n).toArray())
                    return mons.map{ m in TensorGenerator(MonomialGenerator(monomial: m), x) }
                }
                return ModuleObject<Result>(basis: gens)
            },
            differential: { i -> ModuleEnd<Result> in
                let d = self.differential(at: i)
                let cache: CacheDictionary<GridDiagram.Generator, FreeModule<GridDiagram.Generator, P>> = .empty
                return ModuleEnd.linearlyExtend { (t: TensorGenerator<MonomialGenerator<_Un>, GridDiagram.Generator>) -> Result in
                    let (m, x) = t.factors
                    let y = cache.useCacheOrSet(key: x) {
                        d.applied(to: .wrap(x))
                    }
                    return SwiftyKnots.splitMonomials(y).mapGenerators { t in
                        TensorGenerator(m * t.factors.0, t.factors.1)
                    }
                }
            }
        )
    }
}

extension TensorGenerator where A == MonomialGenerator<_Un>, B == GridDiagram.Generator {
    public var algebraicDegree: Int {
        let m = factors.0 // MonomialGenerator<_Un>
        return -(m.monomialDegree.indices.contains(0) ? m.monomialDegree[0] : 0)
    }
    
    public var AlexanderDegree: Int {
        let m = factors.0 // MonomialGenerator<_Un>
        let x = factors.1 // GridDiagram.Generator
        return _Un.totalDegree(exponents: m.monomialDegree) / 2 + x.AlexanderDegree
    }
}

extension GridDiagram {
    public var knotGenus: Int {
        let H = GridComplex.tilde(self).asBigraded { x in x.AlexanderDegree }.homology
        let M = generators.map{ $0.MaslovDegree }.range!
        let A = generators.map{ $0.AlexanderDegree }.range!
        
        for (j, i) in A.reversed() * M.reversed() {
            if !H[i, j].isZero {
                return j
            }
        }
        
        fatalError()
    }
}
