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
import Regex

let dir = "/Users/taketo1024/Projects/swm-kr/data/"
let app = App(storageDir: dir)

app.printResult("3_1", format: .tex)
//app.assertResults("homfly-10")
//app.assertResults("homfly-11")

final class App {
    typealias R = RationalNumber
    typealias Structure = [[Int] : Int]
    
    let storage: Storage
    let maxCrossings: Int
    
    init(storageDir: String, maxCrossings: Int = 11) {
        self.storage = Storage(dir: dir)
        self.maxCrossings = maxCrossings
    }
    
    func load(_ name: String) -> Structure? {
        let file = "\(name).json"
        if !storage.exists(name: file) {
            return nil
        }
        return try? storage.loadJSON(name: file)
    }

    func save(_ name: String, _ str: Structure) {
        let file = "\(name).json"
        try! storage.saveJSON(name: file, object: str)
    }

    func computeAll(_ input: String, skipExisting: Bool = true) {
        let file = "\(dir)/\(input).csv"
        guard let csv = try? CSV(url: URL(fileURLWithPath: file)) else {
            fatalError("Couldn't load \(file)")
        }
        try! csv.enumerateAsDict() { row in
            let target = row["name"]!
            let braidCode = try! JSONDecoder().decode([Int].self, from: row["braid"]!.data(using: .utf8)!)

            let file = "\(target).json"
            if skipExisting && self.storage.exists(name: file) {
                return
            }

            self.compute(target, braidCode)
        }
    }

    func compute(_ target: String, _ braidCode: [Int]) {
        let strands = braidCode.map{ $0.abs }.max()! + 1
        let braid = Braid<anySize>(strands: strands, code: braidCode)
        let K = braid.closure
        
        if K.crossingNumber > self.maxCrossings {
            return
        }
        
        print("target: \(target), n = \(K.crossingNumber)")
        
        let H = KRHomology<R>(K)
        let str = H.structure()
        let raw = str.mapPairs { (g, V) in
            ([g[0], g[1], g[2]], V.rank)
        }
        
        self.save(target, raw)
        
        print(table(raw), "\n")
    }
    
    enum ResultFormat {
        case polynomial, tex, table
    }
    
    func printResults(_ input: String, format: ResultFormat = .polynomial) {
        let file = "\(dir)/\(input).csv"
        guard let csv = try? CSV(url: URL(fileURLWithPath: file)) else {
            fatalError("Couldn't load \(file)")
        }
        try! csv.enumerateAsDict() { row in
            let target = row["name"]!
            self.printResult(target, format: format)
        }
    }
    
    func printResult(_ target: String, format: ResultFormat = .polynomial) {
        guard let str = self.load(target) else {
            return
        }
        
        switch format {
        case .polynomial:
            print(target, ":", asQatPolynomial(str))
        case .tex:
            print(texify(target, asQatPolynomial(str)))
        case .table:
            print(target, "\n", table(str), "\n")
        }
    }

    func assertResults(_ input: String) {
        let file = "\(dir)/\(input).csv"
        guard let csv = try? CSV(url: URL(fileURLWithPath: file)) else {
            fatalError("Couldn't load \(file)")
        }
        try! csv.enumerateAsDict() { row in
            let target = row["name"]!
            let answer = row["HOMFLY"]!
            
            guard let str = self.load(target) else {
                return
            }
            
            let p = self.parseHOMFLY(answer)
            let q = self.asQaPolynomial(str)
            
            if p == q {
                print(target, " : OK")
            } else {
                print(target, " : NG")
                print("\tanswer:", answer)
                print("\tconverted:", p)
                print("\tours     :", q)
            }
        }
    }
    
    private func table(_ structure: Structure) -> String {
        typealias qPoly = KR.qPolynomial<Int>
        let table = structure
            .group { (ijk, _) in
                MultiIndex<_2>(ijk[1], ijk[2])
            }
            .mapValues{ list in
                qPoly(elements: list.map{ (ijk, r) in
                    (ijk[0], r)
                })
            }
        
        return Format.table(
            rows: table.keys.map{ $0[1] }.uniqued().sorted().reversed(), // k
            cols: table.keys.map{ $0[0] }.uniqued().sorted(), // j
            symbol: "k\\j",
            printHeaders: true
        ) { (k, j) -> String in
            let q = table[[j, k]] ?? .zero
            return !q.isZero ? "\(q)" : ""
        }
    }

    private func asQatPolynomial(_ structure: Structure) -> KR.qatPolynomial<Int> {
        .init(elements: structure.map { (g, r) in
            return (MultiIndex(g), r)
        })
    }
    
    public func asQaPolynomial(_ structure: Structure) -> KR.qaPolynomial<Int> {
        .init(elements: structure.map { (g, r) in
            let (i, j, k) = (g[0], g[1], g[2])
            return ([i, j], (-1).pow( (k - j) / 2) * r )
        })
    }
    
    private func texify(_ target: String, _ p: KR.qatPolynomial<Int>) -> String {
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
    
    private func parseHOMFLY(_ poly: String) -> KR.qaPolynomial<Int> {
        typealias P = KR.qaPolynomial<Int>
        
        let z = P(elements: [[1, 0]: 1, [-1, 0]: -1]) // q - q^{-1}
        
        let p0 = #"\^\((-[0-9]+)\)"#
        let r0 = try! Regex(pattern: p0, groupNames: "neg-exp")
        
        let str = r0.replaceAll(in: poly) { m0 in
            let s = m0.group(at: 1)!
            return "^\(s)"
        }
        
        let p1 = #"(\([0-9+\-v\^\*]+\))(\*z(\^(-?[0-9]+))?)?"#
        let p2 = #"([+\-])?(([0-9]+\*)?v(\^-?[0-9]+)?|[0-9]+)"#
        let r1 = try! Regex(pattern: p1, groupNames: "z-terms")
        let r2 = try! Regex(pattern: p2, groupNames: "v-terms")

        return r1.findAll(in: str).sum { m1 -> P in
//            print(m1.matched, m1.subgroups)
            
            let vpoly = m1.subgroups[0]!
            let apoly = r2.findAll(in: vpoly).sum { m2 -> P in
//                print("\t", m2.matched, m2.subgroups)
                let e = m2.subgroups[0].map{ $0 == "-" ? -1 : 1 } ?? 1
                if m2.matched.contains("v") {
                    let c = m2.subgroups[2].flatMap{ Int($0[0 ..< $0.count - 1])} ?? 1
                    let n = m2.subgroups[3].flatMap{ Int($0[1 ..< $0.count]) } ?? 1
                    return P(elements: [[0, n] : e * c])
                } else {
                    let c = m2.subgroups[1].flatMap{ Int($0) }!
                    return P(elements: [[0, 0] : e * c])
                }
            }
            
            let zdeg = m1.subgroups[3].flatMap{ Int($0) } ?? 0
            let qpoly = z.pow(zdeg)

            return apoly * qpoly
        }
    }
}
