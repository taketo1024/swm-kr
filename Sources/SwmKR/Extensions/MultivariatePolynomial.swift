//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/06/16.
//

import SwmCore

extension MultivariatePolynomialType {
    var isLinear: Bool {
        elements.allSatisfy { (e, _) in
            e.total <= 1
        }
    }
    
    var indexOfIndeterminate: Int {
        leadExponent.indices.enumerated().first{ (i, e) in e > 0 }!.offset
    }
    
    func substitute(_ table: [Int : Self]) -> Self {
        let p = RingHom.mapping { i in table[i] }
        return p(self)
    }
}

// TODO: consider AlgebraHom
extension RingHom where Domain: MultivariatePolynomialType, Domain == Codomain {
    public static func mapping(_ f: @escaping (Int) -> Codomain?) -> Self {
        .init { (p: Domain) in
            p.elements.sum { (e, a) in
                a * e.indices.enumerated().multiply { (i, e_i) -> Codomain in
                    if let y = f(i) {
                        return y.pow(e_i) // map x_i^e_i ↦ y^e_i.
                    } else {
                        return .indeterminate(i).pow(e_i) // keep x_i^e_i.
                    }
                }
            }.reduced
        }
    }
}
