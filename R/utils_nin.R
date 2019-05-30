#' @title Not in
#'
#' @name %nin%
#' @rdname not_in
#' @keywords internal
#' @export
#' @importFrom purrr compose
#'
`%nin%` <- compose(`!`, `%in%`)
