// swift-tools-version:5.0
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "SwiftyKnots",
    products: [
        .library(
            name: "SwiftyKnots",
            targets: ["SwiftyKnots"]),
    ],
    dependencies: [
        .package(url: "https://github.com/taketo1024/SwiftyMath.git", .exact("1.0.22")),
        .package(url: "https://github.com/taketo1024/SwiftyMath-homology.git", .exact("1.0.16")),
    ],
    targets: [
        .target(
            name: "SwiftyKnots",
            dependencies: ["SwiftyMath", "SwiftyHomology"],
			path: "Sources/SwiftyKnots"),
        .testTarget(
            name: "SwiftyKnotsTests",
            dependencies: ["SwiftyKnots"]),
        .target(
            name: "SwiftyKnots-Sample",
            dependencies: ["SwiftyMath", "SwiftyHomology", "SwiftyKnots"],
			path: "Sources/Sample"),
    ]
)
