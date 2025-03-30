.make_hex_spots <- function(cdata, fill, coord.multiplier = list(x = 1, y = 1)) {
  assert_that(all(c("x", "y") %in% names(coord.multiplier)))
  
  ## R = circumradius, distance from center to vertex
  ## r = inradius, distance from center to edge midpoint
  r <- 1/2
  R <- (2 / sqrt(3)) * r
  
  spot_positions <- .select_spot_positions(cdata, fill=fill)
  spot_positions <- .adjust_hex_centers(spot_positions)
  
  ## vertices of each hex (with respect to center coordinates)
  ## start at top center, loop clockwise
  vertex_offsets <- data.frame(x.offset=c(0, r, r, 0, -r, -r),
                               y.offset=c(-R, -R/2, R/2, R, R/2, -R/2))
  
  spot_vertices <- .make_spot_vertices(spot_positions, vertex_offsets)
  
  ## Flip to match image orientation
  spot_vertices$x.vertex <- spot_vertices$x.vertex * coord.multiplier$x
  spot_vertices$y.vertex <- spot_vertices$y.vertex * coord.multiplier$y
  
  spot_vertices
}

.adjust_hex_centers <- function(spot_positions) {
  ## R = circumradius, distance from center to vertex
  ## r = inradius, distance from center to edge midpoint
  r <- 1/2
  R <- (2 / sqrt(3)) * r
  
  ## Start at (1-indexed origin)
  spot_positions$x.pos <- spot_positions$x.pos - min(spot_positions$x.pos) + 1
  spot_positions$y.pos <- spot_positions$y.pos - min(spot_positions$y.pos) + 1
  
  ## Shift centers up so rows are adjacent
  spot_positions$y.pos <- spot_positions$y.pos * R * (3/2)
  
  ## Spot columns are offset by row
  ## (i.e. odd rows have odd numbered columns, even rows have even)
  ## Shift centers to the left so columns are adjacent (but hexes stay offset)
  spot_positions$x.pos <- (spot_positions$x.pos + 1) / 2
  
  spot_positions
}

.make_square_spots <- function(cdata, fill="spatial.cluster", scale.factor = 1, offset = 0, coord.multiplier = list(x = 1, y = 1)) {
  assert_that(all(c("x", "y") %in% names(coord.multiplier)))
  
  spot_positions <- .select_spot_positions(cdata, fill=fill)
  
  vertex_offsets <- data.frame(x.offset=c(0, 1, 1, 0),
                               y.offset=c(0, 0, 1, 1))
  vertex_offsets <- vertex_offsets * scale.factor + offset
  
  vertices <- .make_spot_vertices(spot_positions, vertex_offsets)
  
  vertices$x.vertex <- vertices$x.vertex * coord.multiplier$x
  vertices$y.vertex <- vertices$y.vertex * coord.multiplier$y
  
  vertices
}

.select_subspot_positions <- function(cdata, x="spot.col", y="spot.row", fill="spatial.cluster") {
  ## Provide either a column name or vector of labels/values
  assert_that(is.vector(fill) || is.character(fill) || is.factor(fill))
  
  if (is.character(fill) && length(fill) == 1) {
    spot_positions <- cdata[, c(x, y, "subspot.idx", fill)]
    colnames(spot_positions) <- c("x.pos", "y.pos", "subspot.idx", "fill")
  } else if (is.vector(fill) || is.factor(fill)) {
    assert_that(nrow(cdata) == length(fill))
    spot_positions <- cdata[, c(x, y, "subspot.idx")]
    colnames(spot_positions) <- c("x.pos", "y.pos", "subspot.idx")    
    spot_positions$fill <- fill
  }
  
  spot_positions$spot <- rownames(spot_positions)
  
  spot_positions
}

.select_spot_positions <- function(cdata, x="array_col", y="array_row", fill="spatial.cluster") {
  ## Provide either a column name or vector of labels/values
  assert_that(is.vector(fill) || is.character(fill) || is.factor(fill))
  
  ## I think this is the best way to check if something is a string
  if (is.character(fill) && length(fill) == 1) {
    spot_positions <- cdata[, c(x, y, fill)]
    colnames(spot_positions) <- c("x.pos", "y.pos", "fill")    
  } else if (is.vector(fill) || is.factor(fill)) {
    assert_that(nrow(cdata) == length(fill))
    spot_positions <- cdata[, c(x, y)]
    colnames(spot_positions) <- c("x.pos", "y.pos")    
    spot_positions$fill <- fill
  }
  spot_positions$spot <- rownames(spot_positions)
  
  spot_positions
}

.make_spot_vertices <- function(spot_positions, vertex_offsets) {
  spot_vertices <- merge(spot_positions, vertex_offsets)
  spot_vertices$x.vertex <- spot_vertices$x.pos + spot_vertices$x.offset
  spot_vertices$y.vertex <- spot_vertices$y.pos + spot_vertices$y.offset
  
  as.data.frame(spot_vertices)
}

