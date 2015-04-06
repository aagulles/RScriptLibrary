trim.strings <-
function(text, side = "both") {
	if (length(side) != 1) side <- side[1]
	if (is.na(match(side, c("both", "left", "right")))) { side <- "both" }
	if (side == "left") { gsub("\ *$", "", text)	}
	else { if (side == "right") { gsub("^\ *", "", text) }
		 else { gsub("^\ *", "", gsub("\ *$", "", text)) }
	}
}

