//
//  File.swift
//
//
//  Created by Taketo Sano on 2021/06/28.
//

import SwmCore
import SwmKnots
import SwmKR

let dir = "/Users/taketo/Projects/swm/swm-kr/data/"
let app = App(storageDir: dir)

app.logLevel = 3
app.saveResult = true
app.saveProgress = true
app.useBigRational = true
app.useMirror = true

app.computeAll("braid-11")
app.assertResults("homfly-11")
