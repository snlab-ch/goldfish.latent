# direct network
eventsIncrement <- data.frame(
  time = cumsum(
    c(1, 5, 3, 4, 2, 1, 3, 4, 5, 1, 3, 4)),
  sender = sprintf("Actor %d",
                   c(1, 3, 2, 2, 5, 1, 3, 3, 4, 2, 5, 1)),
  receiver = sprintf("Actor %d",
                     c(2, 2, 3, 3, 1, 5, 4, 4, 2, 3, 2, 2)),
  increment =
    c(1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1),
  stringsAsFactors = FALSE
)

actorsEx <- data.frame(
  label = sprintf("Actor %d", 1:5),
  present = c(rep(TRUE, 4), FALSE),
  attr1 = c(9.9, 0.1, 0.5, 0.45, 0.25),
  stringsAsFactors = FALSE
)

compChange <- data.frame(
  node = sprintf("Actor %d", c(5, 4, 4, 1, 5, 1, 5)),
  time = c(10, 12, 17, 26, 26, 30, 30),
  replace = c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE)
)

networkState <- matrix(
  c(0, 3, 0, 0, 0,
    1, 0, 1, 1, 0,
    0, 0, 0, 1, 0,
    0, 0, 1, 0, 0,
    0, 0, 0, 0, 0),
  nrow = 5, ncol = 5, byrow = TRUE,
  dimnames = list(sprintf("Actor %d", 1:5),
                  sprintf("Actor %d", 1:5))
)

# exogenous network
eventsExogenous <- data.frame(
  time =
    c(7, 14, 15, 18, 18, 25, 25),
  sender = sprintf("Actor %d",
                   c(4,  2,  5,  4,  4,  1,  3)),
  receiver = sprintf("Actor %d",
                     c(2,  3,  1,  5,  2,  3,  5)),
  increment =
    c(1,  1,  3,  1, -1,   2, 3),
  stringsAsFactors = FALSE
)

networkExog <- matrix(
  c(0, 0, 0, 1, 0,
    0, 0, 0, 0, 0,
    0, 2, 0, 0, 0,
    1, 0, 0, 0, 0,
    1, 2, 0, 0, 0),
  nrow = 5, ncol = 5, byrow = TRUE,
  dimnames = list(sprintf("Actor %d", 1:5),
                  sprintf("Actor %d", 1:5))
)

# defining objects
actorsEx <- goldfish::defineNodes(actorsEx) |>
  goldfish::linkEvents(changeEvent = compChange, attribute = "present")

networkState <- goldfish::defineNetwork(
  matrix = networkState, nodes = actorsEx,
  directed = TRUE) |>
  goldfish::linkEvents(changeEvent = eventsIncrement, nodes = actorsEx)
depNetwork <- goldfish::defineDependentEvents(
  events = eventsIncrement,
  nodes = actorsEx,
  defaultNetwork = networkState)

# define goldfish objects
networkExog <- goldfish::defineNetwork(
  matrix = networkExog,
  nodes = actorsEx, directed = TRUE) |>
  goldfish::linkEvents(changeEvent = eventsExogenous, nodes = actorsEx)

