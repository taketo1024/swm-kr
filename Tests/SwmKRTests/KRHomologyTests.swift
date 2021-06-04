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
    
    override func setUp() {}
    
    override func tearDown() {}
    
    func testTrefoil() {
        typealias R = 𝐐
        
        let K = Link.load("3_1")!
        let H = KRHomology<R>(K)
        let str = H.structure()
        
        XCTAssertEqual(str.count, 3)
        XCTAssertEqual(str[[0, -4, 2]]!.rank, 1)
        XCTAssertEqual(str[[-2, -2, 2]]!.rank, 1)
        XCTAssertEqual(str[[2, -2, -2]]!.rank, 1)
    }
    
    func testFigureEight() {
        typealias R = 𝐐
        
        let K = Link.load("4_1")!
        let H = KRHomology<R>(K)
        let str = H.structure()
        
        XCTAssertEqual(str.count, 5)
        XCTAssertEqual(str[[2, 0, -2]]!.rank, 1)
        XCTAssertEqual(str[[-2, 0, 2]]!.rank, 1)
        XCTAssertEqual(str[[0, 2, -2]]!.rank, 1)
        XCTAssertEqual(str[[0, -2, 2]]!.rank, 1)
        XCTAssertEqual(str[[0, 0, 0]]!.rank, 1)
    }
}
