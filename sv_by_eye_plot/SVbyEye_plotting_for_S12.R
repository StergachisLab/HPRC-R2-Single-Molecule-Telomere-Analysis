

# CORRECT: Just point directly to your file
paf.file <- "telo_100kb.ava.alignment.for_S12.paf"

# Then run readPaf
paf.table <- readPaf(
  paf.file = paf.file,
  include.paf.tags = TRUE, 
  restrict.paf.tags = "cg"
)
paf.table
unique(paf.table$t.name)
flip_id <- "HG03098#1#CM092238.1:93010326-93217697"
#---- flip when it appears as query ----
iq <- paf.table[["q.name"]] == flip_id
qs <- paf.table[["q.start"]][iq]; qe <- paf.table[["q.end"]][iq]; ql <- paf.table[["q.len"]][iq]
paf.table[["q.start"]][iq] <- ql - qe
paf.table[["q.end"]][iq]   <- ql - qs

# ---- flip when it appears as target ----
it <- paf.table[["t.name"]] == flip_id
ts <- paf.table[["t.start"]][it]; te <- paf.table[["t.end"]][it]; tl <- paf.table[["t.len"]][it]
paf.table[["t.start"]][it] <- tl - te
paf.table[["t.end"]][it]   <- tl - ts

#  flip strand too so direction coloring makes sense
paf.table[["strand"]][iq | it] <- ifelse(paf.table[["strand"]][iq | it] == "+", "-", "+")

paf.table <- paf.table[(paf.table[["t.end"]] - paf.table[["t.start"]]) > 10000, ]
seqnames.order <- c("NA18570#2#CM087682.1:1-109168","HG03050#1#CM098768.1:1-108248",
                    "HG03139#2#CM091842.1:1-112107", "HG03225#1#CM089305.1:1-105605", 
                    "NA19835#1#CM094029.1:1-107376","HG03098#1#CM092238.1:93110326-93217697",
                    "HG00099#1#CM087326.1:1-112798", "HG00099#2#CM087370.1:1-111965", 
                    "HG00140#1#CM087120.1:1-111680","HG00253#1#CM094466.1:1-104387", 
                    "HG00272#2#CM094225.1:1-109223","HG02392#2#CM099896.1:1-111146",
                    "NA18505#2#CM089876.1:1-107285")

length(seqnames.order)
length(unique(paf.table$t.name))

plotAVA(
  paf.table = paf.table,
  color.by = "direction",
  seqnames.order = seqnames.order
)
