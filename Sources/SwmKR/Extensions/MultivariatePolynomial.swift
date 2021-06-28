//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/06/16.
//

import SwmCore

extension MultiIndex {
    var firstNonZeroIndex: Int? {
        indices.enumerated().first{ (i, e) in e > 0 }?.offset
    }
}

extension MultivariatePolynomialType {
    var indeterminateOfPrimaryExclusion: Int? {
        if !elements.allSatisfy({ (e, _) in
            e.total <= 1
        }) {
            return nil
        }
        return elements.first{ (e, a) in !a.isZero }?.key.firstNonZeroIndex
    }
    
    var indeterminateOfSecondaryExclusion: Int? {
        if !elements.allSatisfy({ (e, _) in
            e.total <= 2
        }) {
            return nil
        }
        return elements.first{ (e, a) in !a.isZero && e.indices.contains(2) }?.key.firstNonZeroIndex
    }
    
    func substitute(_ table: [Int : Self]) -> Self {
        if table.isEmpty {
            return self
        } else {
            let p = RingHom.mapping { i in table[i] }
            return p(self)
        }
    }
    
    private func highestExponent(as i: Int) -> Exponent {
        elements
            .exclude{ (e, a) in a.isZero || e[i] == 0 }
            .map{ (e, _) in e }
            .max{ e in e[i] }
        ?? .zero
    }
    
    func divide(by g: Self, as i: Int) -> (quotient: Self, remainder: Self) {
        let e0 = g.highestExponent(as: i)
        let a0 = g.coeff(e0)
        
        var q = Self.zero
        var r = self
        
        while !r.isZero {
            let e1 = r.highestExponent(as: i)
            let e = e1 - e0
            if e.indices.contains(where: {$0 < 0} ) {
                break
            }
            
            let a = r.coeff(e1) * a0.inverse!
            let x = a * Self.monomial(withExponents: e)
            
            r = r - x * g
            q = q + x
        }
        
        return (q, r)
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
                        return .indeterminate(i, exponent: e_i) // keep x_i^e_i.
                    }
                }
            }
        }
    }
}
