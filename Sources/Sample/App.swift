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
        case polynomial, table, texPolynomial, texTable
    }
    
    func printResults(_ input: String, format: ResultFormat = .table) {
        let file = "\(dir)/\(input).csv"
        guard let csv = try? CSV(url: URL(fileURLWithPath: file)) else {
            fatalError("Couldn't load \(file)")
        }
        try! csv.enumerateAsDict() { row in
            let target = row["name"]!
            self.printResult(target, format: format)
        }
    }
    
    func printResult(_ target: String, mirror: Bool = false, format: ResultFormat = .table) {
        guard let str = self.load(target) else {
            return
        }
        
        func _printResult(_ target: String, _ str: Structure) {
            switch format {
            case .table:
                print("\(target)\n\(table(str))\n")
            case .polynomial:
                print("\(target) : \(PoincarePolynomial(str))")
            case .texPolynomial:
                print(texPolynomial(target, str))
            case .texTable:
                print(texTable(target, str))
            }
        }
        
        mirror ? _printResult("m(\(target))", self.mirror(str)) : _printResult(target, str)
    }
    
    func isDistinctPair(_ K1: String, _ K2: String) -> Bool {
        let str1 = load(K1)!
        let str2 = load(K2)!
        
        return !(str1 == str2 || str1 == mirror(str2))
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
            let q = self.HOMFLYPolynomial(str)
            
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
    
    func findKRThick(_ input: String) -> [String] {
        let file = "\(dir)/\(input).csv"
        guard let csv = try? CSV(url: URL(fileURLWithPath: file)) else {
            fatalError("Couldn't load \(file)")
        }
        
        var res: [String] = []
        try! csv.enumerateAsDict() { row in
            let target = row["name"]!
            guard let str = self.load(target) else {
                return
            }
            if !self.isKRThin(str) {
                res.append(target)
            }
        }
        
        return res
    }
    
    private func mirror(_ structure: Structure) -> Structure {
        structure.mapPairs{ (g, r) in ([-g[0], -g[1], -g[2]], r)}
    }
    
    private func isKRThin(_ structure: Structure) -> Bool {
        let s = structure.keys.map{ g in g[0] + g[1] + g[2]}
        return Set(s).count == 1
    }
    
    public func HOMFLYPolynomial(_ structure: Structure) -> KR.qaPolynomial<Int> {
        .init(elements: structure.map { (g, r) in
            let (i, j, k) = (g[0], g[1], g[2])
            return ([i, j], (-1).pow( (k - j) / 2) * r )
        })
    }
    
    public func PoincarePolynomial(_ structure: Structure) -> KR.tqaPolynomial<Int> {
        .init(elements: structure.map { (g, r) in
            let (i, j, k) = (g[0], g[1], g[2])
            return ([(k - j) / 2, i, j], r)
        })
    }
    
    private func tableData(_ structure: Structure) -> ([MultiIndex<_2> : KR.qPolynomial<Int>], [Int], [Int]) {
        typealias qPoly = KR.qPolynomial<Int>
        let table = structure
            .group { (g, _) in
                MultiIndex<_2>(g[1], g[2])
            }
            .mapValues{ list in
                qPoly(elements: list.map{ (g, r) in
                    (g[0], r)
                })
            }
        
        let js = table.keys.map{ $0[0] }.uniqued().sorted()
        let ks = table.keys.map{ $0[1] }.uniqued().sorted()
        return (table, js, ks)
    }
    
    private func table(_ structure: Structure) -> String {
        let (table, js, ks) = tableData(structure)
        return Format.table(
            rows: ks.reversed(),
            cols: js,
            symbol: "k\\j",
            printHeaders: true
        ) { (k, j) -> String in
            let q = table[[j, k]] ?? .zero
            return !q.isZero ? "\(q)" : ""
        }
    }

    private func texPolynomial(_ target: String, _ structure: Structure) -> String {
        return structure.mapPairs { (g, r) -> ([Int], Int) in
            let (i, j, k) = (g[0], g[1], g[2])
            let h = (k - j) / 2
            return ([h, i, j], r)
        }.sorted { (f, g) in
            f.0.reversed().lexicographicallyPrecedes(g.0.reversed()) // a > q > t
        }.map { (g, r) -> String in
            let mon = zip(["t", "q", "a"], g).map { (x, i) in
                i == 0 ? "" : "\(x)^{\(i)}"
            }.joined()
            return mon == "" ? "\(r)" : r == 1 ? mon : "\(r)\(mon)"
        }.joined(separator: " + ")
    }
    
    private func texTable(_ target: String, _ structure: Structure) -> String {
        let (table, js, ks) = tableData(structure)
        let body = ks.reversed().map { k -> String in
            "$\(k)$ & " + js.map { j in
                let qpoly = table[[j, k]] ?? .zero
                return qpoly.terms.sorted{ $0.leadExponent }.map { term -> String in
                    let n = term.leadCoeff
                    let i = term.leadExponent
                    let x = i == 0 ? "\(n)" : n == 1 ? "q^{\(i)}" : "\(n)q^{\(i)}"
                    return "$\(x)$"
                }.joined(separator: " + ")
            }.joined(separator: " & ")
        }.joined(separator: " \\\\\n") + " \\\\"
        
        return """
\\begin{minipage}{\\textwidth}
\\item $\(target)$ \\vspace{0.5em} \\\\
\\begin{tabular}{l|\(Array(repeating: "l", count: js.count).joined())}
$k \\setminus j$ & \(js.map{ j in "$\(j)$" }.joined(separator: " & ")) \\\\
\\hline
\(body)
\\end{tabular}
\\vspace{1em}
\\end{minipage}
%
"""
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
