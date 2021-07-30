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

let dir = "/Users/taketo1024/Projects/swm-kr/data/"
let app = App(storageDir: dir)

app.printResults("braid-11", format: .table)
