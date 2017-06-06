import Foundation

public protocol Module: AdditiveGroup {
    associatedtype R: Ring
    static func * (r: R, m: Self) -> Self
    static func * (m: Self, r: R) -> Self
}

public protocol Submodule: Module, AdditiveSubgroup {
    associatedtype Super: Module
}

public extension Submodule where R == Super.R {
    static func * (r: R, a: Self) -> Self {
        return Self.init(r * a.asSuper)
    }
    
    static func * (a: Self, r: R) -> Self {
        return Self.init(a.asSuper * r)
    }
}

public protocol ProductModuleType: Module, AdditiveProductGroup {
    associatedtype Left: Module
    associatedtype Right: Module
}

public extension ProductModuleType where Left.R == R, Right.R == R {
    static func * (r: R, a: Self) -> Self {
        return Self.init(r * a._1, r * a._2)
    }
    
    static func * (a: Self, r: R) -> Self {
        return Self.init(a._1 * r, a._2 * r)
    }
}

public struct ProductModule<M1: Module, M2: Module>: ProductModuleType where M1.R == M2.R {
    public typealias Left = M1
    public typealias Right = M2
    public typealias R = M1.R
    
    public let _1: M1
    public let _2: M2
    
    public init(_ m1: M1, _ m2: M2) {
        self._1 = m1
        self._2 = m2
    }
}

public protocol QuotientModuleType: Module, AdditiveQuotientGroup {
    associatedtype Sub: Submodule
}

public extension QuotientModuleType where R == Sub.R, R == Sub.Super.R {
    public static func isEquivalent(_ a: Base, _ b: Base) -> Bool {
        return Sub.contains( a - b )
    }
    
    static func * (r: R, a: Self) -> Self {
        return Self.init(r * a.representative)
    }
    
    static func * (a: Self, r: R) -> Self {
        return Self.init(a.representative * r)
    }
    
    public var hashValue: Int {
        return representative.hashValue // must assure `representative` is unique.
    }
}

public struct QuotientModule<M: Module, S: Submodule>: QuotientModuleType where M == S.Super, M.R == S.R {
    public typealias R = M.R
    public typealias Sub = S
    
    internal let m: M
    
    public init(_ m: M) {
        self.m = m // TODO reduce
    }
    
    public var representative: M {
        return m
    }
}
