//
//  File.swift
//  
//
//  Created by Taketo Sano on 2021/08/16.
//

import SwmCore
import SwmMatrixTools
import SwmHomology
import SwmxBigInt

extension BigRational: ComputationalField {
    public typealias ComputationalMatrixImpl = DefaultSparseMatrixImpl<Self> // TODO
    public typealias ComputationalSparseMatrixImpl = DefaultSparseMatrixImpl<Self>
    
    public var computationalWeight: Double {
        isZero ? 0 : Double(max(numerator.abs, denominator))
    }
}

extension BigRational: HomologyCalculatable{}

