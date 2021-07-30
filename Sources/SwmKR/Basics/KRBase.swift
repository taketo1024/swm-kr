//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/05/31.
//

import SwmCore
import SwmHomology
import SwmKhovanov
import SwmKnots

public struct KR {
    public struct Indeterminates {
        public struct x: PolynomialIndeterminate {
            public static let degree = 2
            public static var symbol = "x"
        }
        public typealias xn = EnumeratedPolynomialIndeterminates<x, anySize>
        
        public struct t: PolynomialIndeterminate {
            public static let symbol = "t"
        }
        public struct q: PolynomialIndeterminate {
            public static let symbol = "q"
        }
        public struct a: PolynomialIndeterminate {
            public static let symbol = "a"
        }
        
        public typealias qa = BivariatePolynomialIndeterminates<q, a>
        public typealias tqa = TrivariatePolynomialIndeterminates<t, q, a>
    }

    public typealias Grading = MultiIndex<_3>
    public typealias EdgeRing<R: Ring> = MultivariatePolynomial<R, Indeterminates.xn>
    public typealias BaseModule<R: Ring> = LinearCombination<R, MonomialAsGenerator<Indeterminates.xn>>
    public typealias HorizontalModule<R: Ring> = IndexedModule<Cube.Coords, BaseModule<R>>
    public typealias TotalModule<R: Ring> = IndexedModule<Cube.Coords, HorizontalModule<R>>
    
    public typealias qPolynomial<R: Ring>  = LaurentPolynomial<R, Indeterminates.q>
    public typealias qaPolynomial<R: Ring> = MultivariateLaurentPolynomial<R, Indeterminates.qa>
    public typealias tqaPolynomial<R: Ring> = MultivariateLaurentPolynomial<R, Indeterminates.tqa>
    
    public typealias Structure<R: Ring> = [KR.Grading : ModuleStructure<KR.TotalModule<R>>]

    static func baseGrading(link L: Link, hCoords: Cube.Coords, vCoords: Cube.Coords) -> KR.Grading {
        (0 ..< L.crossingNumber).sum { i -> KR.Grading in
            switch (L.crossings[i].crossingSign, hCoords[i], vCoords[i]) {
            case (+1, 0, 0):
                return [2, -2, -2]
            case (+1, 1, 0):
                return [0, 0, -2]
            case (+1, 0, 1):
                return [0, -2, 0]
            case (+1, 1, 1):
                return [0, 0, 0]
                
            case (-1, 0, 0):
                return [0, -2, 0]
            case (-1, 1, 0):
                return [0, 0, 0]
            case (-1, 0, 1):
                return [0, -2, 2]
            case (-1, 1, 1):
                return [-2, 0, 2]
                
            default:
                fatalError("impossible")
            }
        }
    }

    struct EdgeConnection<R: Ring>: CustomStringConvertible {
        let ik: EdgeRing<R>
        let il: EdgeRing<R>
        
        var description: String {
            "\((ik: ik, il: il))"
        }
    }
}
