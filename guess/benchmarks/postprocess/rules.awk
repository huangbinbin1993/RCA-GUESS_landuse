function o(x) {
  print $1, $2, x
}

BEGIN {
  pft = 20
  lat = 2
  trees = 19
  grass = 18
  tr = 16
  te = 17
  bt = 15
  total = 14
  trbe = 9
  trbr = 11
}

{
if ($trees > .5 && $pft == "BNS" && $bt > .8*$trees) o(1)
else if ($trees > .5 && ($pft == "BNE" || $pft == "IBS") && $bt > .8*$trees) o(2)
else if ($trees > 2.5 && $trees*.2 < $bt && $bt < .8*$trees && .2*$trees < $te && $te < .8*$trees) o(3)
else if ($pft == "TeBS" && $trees > 2.5 && $te > .8*$trees) o(5)
else if ($pft == "TeBE" && $trees > 2.5 && $te > .8*$trees) o(6)
else if ($trees > 2.5 && $te > .8*$trees) o(7)
else if ($trees > 2.5 && $tr > .5*$trees && $trbe <=.6*$trees && $trbr <= .6*$trees) o(8)
else if ($trees > 2.5 && $trbe > .6*$trees && $pft == "TrBE") o(9)
else if ($trees > 2.5 && $trbr > .6*$trees && $pft == "TrBR") o(10)
else if ($trees < .5 && $total > .2 && $lat > 54 && $trees < .5*$grass) o(18)
else if ($trees < 2.5 && $trees > .5 && $grass < $trees) o(15)
else if ($total > 3 && $trees > .5 && $trees < 2.5) o(11)
else if ($total <= 3 && $trees > .5 && $trees < 2.5) o(12)
else if ($grass > 3 && $trees < .5) o(13)
else if ($grass > .5 && $trees < .2) o(14)
else if ($trees < .5 && $total > .2) o(16)
else if ($total <= .2) o(17)
}
