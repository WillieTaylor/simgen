  BEGIN{RS="@"}
  /Taylor, W\. R\./ {print "@"$0}
