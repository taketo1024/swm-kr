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
    
    private func highestExponent(as i: Int) -> Exponent {
        elements
            .exclude{ (e, a) in a.isZero || e[i] == 0 }
            .map{ (e, _) in e }
            .max{ e in e[i] }
        ?? .zero
    }
    
    func divide(by g: Self, as i: Int) -> Self {
//        print("\t\tdivide \(self) / \(g)")

        var f = self
        var q = Self.zero
        let e0 = g.highestExponent(as: i)
        let a0 = g.coeff(e0)
        
        assert(!e0.isZero)
        
        while !f.isZero {
//            print("\t\t\tdivide \(f) / \(g)")
            let e1 = f.highestExponent(as: i)
            let e = e1 - e0
            assert(e1 != .zero)
            assert(e.indices.allSatisfy{ $0 >= 0 })
            
            let a = f.coeff(e1) * a0.inverse!
            let x = a * Self.monomial(withExponents: e)
            
            f = f - x * g
            q = q + x
        }
        
        assert(self == q * g)
        return q
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
