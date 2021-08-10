//
//  File.swift
//
//
//  Created by Taketo Sano on 2021/06/28.
//

import Foundation
import SwmCore
import SwmKnots
import SwmKR
import Fortify

let dir = "/Users/taketo/Projects/swm/swm-kr/data/"
let app = App(storageDir: dir, maxCrossings: 12)

do {
    try Fortify.protect {
        app.computeAll("braid-11")
    }
} catch {
    print("error: \(error)")
}
