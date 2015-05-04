"use strict"

module.exports = splatPoints

var dirichlet = require("dirichlet")
var ndarray = require("ndarray")
var doSplat = require("cwise")({
  args: [
    "array",
    { "offset" : [0, 1], "array": 0 },
    "array",
    "scalar",
    "scalar",
    "scalar"],
  pre: function(px, py, w, out, radius, dirichlet) {
    var os = out.shape
    this.nx = os[0]|0
    this.ny = os[1]|0
    this.ir = Math.ceil(radius)|0
    this.min = Math.min
    this.max = Math.max
    this.floor = Math.floor
  },
  body: function(x, y, w, out, radius, dirichlet) {
    var ix = this.floor(x)|0
    var iy = this.floor(y)|0
    var x0 = this.min(this.max(ix - this.ir, 0), this.nx)
    var x1 = this.min(this.max(ix + this.ir, 0), this.nx)
    var y0 = this.min(this.max(iy - this.ir, 0), this.ny)
    var y1 = this.min(this.max(iy + this.ir, 0), this.ny)

    //Splat point
    for(var i=x0; i<x1; ++i) {
      var wx = w * dirichlet(this.nx, x-i)
      for(var j=y0; j<y1; ++j) {
        var wy = wx * dirichlet(this.ny, y-j)
        out.set(i,j, out.get(i,j) + wy)
      }
    }
  }
})

//Splat a list of points to a grid with the given weights
function splatPoints(out, points, weights, radius) {
  doSplat(
    points.hi(points.shape[0], 1),
    ndarray(weights.data,
      [weights.shape[0], 1],
      [weights.stride[0], 0],
      weights.offset),
    out,
    radius,
    dirichlet)
}
