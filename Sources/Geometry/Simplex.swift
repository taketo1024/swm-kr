//
//  Simplex.swift
//  SwiftyAlgebra
//
//  Created by Taketo Sano on 2017/05/03.
//  Copyright © 2017年 Taketo Sano. All rights reserved.
//

import Foundation

// MEMO: 'un'ordered set of vertices

public struct Simplex: FreeModuleBase, CustomStringConvertible {
    public let vertices: [Vertex]
    private let verticesSet: Set<Vertex>
    private let id: String
    
    public var dim: Int {
        return vertices.count - 1
    }
    
    internal init<S: Sequence>(_ vertices: S) where S.Iterator.Element == Vertex {
        self.vertices = vertices.sorted().unique()
        self.verticesSet = Set(self.vertices)
        
        self.id = "(\(self.vertices.map{$0.description}.joined(separator: ", ")))"
    }
    
    public func face(_ index: Int) -> Simplex {
        let vs = (0 ... dim).filter({$0 != index}).map{vertices[$0]}
        return Simplex(vs)
    }
    
    public func faces() -> [Simplex] {
        if dim == 0 {
            return []
        } else {
            return (0 ... dim).map{ face($0) }
        }
    }
    
    public func contains(_ s: Simplex) -> Bool {
        return s.verticesSet.isSubset(of: self.verticesSet)
    }
    
    public func allSubsimplices() -> [Simplex] {
        var queue = [self]
        var i = 0
        while(i < queue.count) {
            let s = queue[i]
            if s.dim > 0 {
                queue += queue[i].faces()
            }
            i += 1
        }
        return queue.unique()
    }
    
    public var hashValue: Int {
        return id.hashValue
    }
    
    public var description: String {
        return id
    }
    
    public static func ==(a: Simplex, b: Simplex) -> Bool {
        return a.id == b.id // should be `a.verticesSet == b.verticesSet` but for performance.
    }
}