module LickometerAnalysis
using DelimitedFiles
using PyPlot
using Peaks
using Images
using Statistics

export x_ms2min, x_ms2sec, setthresh, roll1d, baselinecorrection, findlicks, filterlicks, findlickratePerBout, min2time, time2min, find_bout_boundary, lickcount

# Write your package code here.
"""convert milliseconds to miniutes"""
function x_ms2min(cycles, waittime)
  xmin = 0:(waittime/1000/60):(cycles-1)*(waittime/1000/60) #Convert to minutes
  return(xmin)
end

"""convert milliseconds to seconds"""
function x_ms2sec(cycles, waittime)
  xsec = 0:(waittime/1000):(cycles-1)*(waittime/1000) #Convert to minutes
  return(xsec)
end

"""threshold is the `frac` of the `maxidx`-th highest value in the `dat` vector"""
function setthresh(dat, maxidx, frac)
  thresh = sort(dat, rev = true)[maxidx]*frac
  return(thresh)
end

function roll1d(f, v, win)
  output = zeros(length(v)-win)
  for (i, val) in enumerate(Int(1-win/2):Int(length(v)-win-win/2))
    output[i] = f(v[val:val+win-1])
  end
  output
end

function baselinecorrection(dat, win, filt)
  tmp = padarray(dat, Pad(:symmetric, win))
  val  = roll1d(filt, tmp, win*2)
  datn = dat .- val
  return(datn)
end

"""Finding local maxima that are higher than the threshold"""
function findlicks(v, thresh, win)
  vout = Vector{Int}()
  v_jitter = v .+ rand(length(v)).*1e-10
  localmax_idx, localmax_val = findmaxima(v_jitter, win) ## Critial factor = window size
  for (i, val) in enumerate(localmax_val)
    if val > thresh
      push!(vout, localmax_idx[i])
    end
  end
  vout
end

"""Finding local maxima that are higher than the threshold with `filterlicks`"""
function findlicks(vec, thresh, win, filterwinsize)
  vout = findlicks(vec, thresh, win)
  lickidx = filterlicks(vec, thresh, vout, filterwinsize)
  return(lickidx)
end

""" Find local maxima whose neighbors are lower than the threshold"""
function filterlicks(vec::AbstractVector, thresh::Number, idx::Number, winsize::Number)
  vec_thresh = vec .> thresh
  return(!all(vec_thresh[idx-winsize:idx+1]))
end

function filterlicks(vec::AbstractVector, thresh::Number, lickidx::AbstractVector, winsize::Number)
  keep = falses(length(lickidx))
  for (n, i) in enumerate(lickidx)
    keep[n] = filterlicks(vec, thresh, i, winsize)
  end
  return(lickidx[keep])
end

function findlickratePerBout(lickidx, waittime, boutinterval_thresh)
  bout_boundary = findall(diff(lickidx) .> boutinterval_thresh)
  dts = [(lickidx[bout_boundary[i+1]] - lickidx[bout_boundary[i]+1])*waittime for i in 1:length(bout_boundary)-1]
  nlicksperbout = [bout_boundary[i+1]-bout_boundary[i]+1+1 for i in 1:length(bout_boundary)-1]
  lickratePerBout = nlicksperbout./dts*1000
  lickratePerBout
end

""" hour:min:sec to minutes. `offset` in second"""
function min2time(minute; offset = 0, scale = 1)
  sec = minute*60
  sec = sec-offset
  minute = sec/60
  minute = minute * scale
  hour = minute/60
  h = floor(Int, hour) 
  m = hour%1*60
  s = m%1*60
  m = floor(Int, m)
  return(h, m, s)
end  

""" minute to hour:min:sec `offset` in second"""
function time2min(h, m, s; offset = 0)
  s = s-offset
  m = m + s/60
  return(h*60+m)
end

"""find lick boundary: A bout is consist of multiple licks"""
function find_bout_boundary(lickidx, boutinterval_thresh)
  d_lickidx = diff(lickidx)
  bout_boundary = findall(d_lickidx .> boutinterval_thresh)
  v = Vector{UnitRange}(undef, length(bout_boundary));
  for i in eachindex(bout_boundary)
    if i == 1
      v[i] = lickidx[i]:lickidx[bout_boundary[i]]
    else
      v[i] = d_lickidx[bout_boundary[i-1]] + lickidx[bout_boundary[i-1]]: lickidx[bout_boundary[i]]
    end
  end
  v
end

"""
Count licks within a give time window.
Time window can be set by `timefirst` and `timelast` in ms, sec, or min, or else depending on the unit of `x_axis`.
"""
function lickcount(x_axis, lickidx, timefirst, timelast)
  keep = x_axis .>= timefirst .&& x_axis .<= timelast #Time window from 12am - 8am.
  firstidx = findfirst(keep .== 1)
  lastidx = findlast(keep .== 1)
  valout = sum(lickidx .>= firstidx .&& lickidx .<= lastidx)
  return(valout)
end

end
