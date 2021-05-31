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
    }

    public typealias Grading = MultiIndex<_3>
    public typealias EdgeRing<R: Ring> = MultivariatePolynomial<R, Indeterminates.xn>
    public typealias BaseModule<R: Ring> = LinearCombination<R, MonomialAsGenerator<Indeterminates.xn>>
    public typealias HorizontalModule<R: Ring> = GradedModule<Cube.Coords, BaseModule<R>>
    public typealias TotalModule<R: Ring> = GradedModule<Cube.Coords, HorizontalModule<R>>

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

