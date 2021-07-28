//
//  File.swift
//
//
//  Created by Taketo Sano on 2021/06/28.
//

import Foundation
import SwmCore
import SwmHomology
import SwmKnots
import SwmKR
import SwiftCSV

typealias R = RationalNumber

let storage = Storage(dir: "/Users/taketo/Projects/swm/swm-kr/data/")

saveData()
//printList()

func saveData() {
    let csv: CSV = try! CSV(url: URL(fileURLWithPath: "/Users/taketo/Projects/swm/swm-kr/data/data.csv"))
    try! csv.enumerateAsDict() { row in
        let target = row["name"]!
        let code = try! JSONDecoder().decode([Int].self, from: row["braid"]!.data(using: .utf8)!)
        let strands = code.map{ $0.abs }.max()! + 1
        let braid = Braid<anySize>(strands: strands, code: code)
        compute(target, braid)
    }

//    for n in 3 ... 10 {
//        for i in 1 ... 200 {
//            let target = "\(n)_\(i)"
//            let file = "\(target).json"
//
//            if storage.exists(name: file) {
//                continue
//            }
//
//            guard let b = Braid.load(target) else {
//                break
//            }
//
//            compute(target, b)
//        }
//    }
}

func compute(_ target: String, _ b: Braid<anySize>) {
    let file = "\(target).json"
    if storage.exists(name: file) {
        return
    }
    
    let K = b.closure
    
    if K.crossingNumber >= 12 {
        return
    }
    
    print(Date())
    print("target:", target, "/ n =", K.crossingNumber)

    let H = KRHomology<R>(K)
    let str = H.structure()
    try! storage.saveJSON(name: "\(target).json", object: str)

    print(H.table(), "\n")
    print(KR.asQatPolynomial(str), "\n")
}

func printList() {
    for n in 3 ... 10 {
        for i in 1 ... 200 {
            let target = "\(n)_\(i)"
            let file = "\(target).json"
            
            if !storage.exists(name: file) {
                continue
            }

            let str: KR.Structure<R> = try! storage.loadJSON(name: file)
            let p = KR.asQatPolynomial(str)
            let tex = texify(target, p)
            print(tex)
        }
    }
}

func texify(_ target: String, _ p: KR.qatPolynomial<Int>) -> String {
    func aqt(_ p: KR.qatPolynomial<Int>) -> [Int] {
        let g = p.leadExponent.indices
        return [g[1], g[0], g[2]]
    }
    let poly = p.terms.sorted { (f, g) in
        aqt(f).lexicographicallyPrecedes(aqt(g))
    }.map { f -> String in
        let e = f.leadExponent.indices
        let (i, j, h) = (e[0], e[1], e[2])
        return "t^{\(h)}q^{\(i)}a^{\(j)}"
    }.joined(separator: " + ")
    return "$\(target)$ & $\(poly)$ \\\\"
}
