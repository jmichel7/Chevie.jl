function famg4()
  g4=crg(4)
  classinfo(g4)
  g4.classinfo.classnames=["1","z","g_4","g_6","g_6^4","g_6^2","g_6^5"]
  drinfeld_double(g4;pivotal=(g4(1,2)^3,[E(3),E(3)]))
end
