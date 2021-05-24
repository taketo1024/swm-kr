//
//  KHTests.swift
//  SwiftyKnots
//
//  Created by Taketo Sano on 2018/04/04.
//

import XCTest
import SwmCore
import SwmKnots
import SwmHomology
@testable import SwmKR

class KRHomologyTests: XCTestCase {
    
    typealias R = 𝐐
    typealias KR = KRHomology<R>
    
    override func setUp() {}
    
    override func tearDown() {}
    
    func testTrefoil() {
        let K = Link.load("3_1")!
        let H = KR(K)
        print(H.structure())
    }
}
