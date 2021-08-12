// swift-tools-version:5.3
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "swm-kr",
    products: [
        .library(
            name: "SwmKR",
            targets: ["SwmKR"]
        ),
    ],
    dependencies: [
        .package(
			url: "https://github.com/taketo1024/swm-core.git",
			from:"1.2.9"
//            path: "../swm-core/"
		),
        .package(
			url: "https://github.com/taketo1024/swm-knots.git",
			from: "1.2.0"
//            path: "../swm-knots/"
		),
        .package(
            url: "https://github.com/taketo1024/swm-homology.git",
            from: "1.3.3"
//            path: "../swm-homology/"
		),
        .package(
            url: "https://github.com/taketo1024/swm-khovanov.git",
            from: "1.1.6"
//            path: "../swm-khovanov/"
        ),
        .package(
            url: "https://github.com/taketo1024/swmx-bigint.git",
            from: "1.0.0"
//            path: "../swmx-bigint/"
        ),
        .package(
            url: "https://github.com/swiftcsv/SwiftCSV.git",
            from: "0.6.0"
        ),
        .package(
            url: "https://github.com/crossroadlabs/Regex.git",
            from: "1.2.0"
        ),
        .package(url: "https://github.com/johnno1962/Fortify", from: "2.1.4")
    ],
    targets: [
        .target(
            name: "SwmKR",
            dependencies: [
                .product(name: "SwmCore", package: "swm-core"),
                .product(name: "SwmKnots", package: "swm-knots"),
                .product(name: "SwmHomology", package: "swm-homology"),
                .product(name: "SwmKhovanov", package: "swm-khovanov"),
            ]
        ),
        .target(
            name: "Sample",
            dependencies: [
                "SwmKR",
                "SwiftCSV",
                "Regex",
                "Fortify",
                .product(name: "SwmxBigInt", package: "swmx-bigint"),
            ]),
        .testTarget(
            name: "SwmKRTests",
            dependencies: ["SwmKR"]
		),
    ]
)
