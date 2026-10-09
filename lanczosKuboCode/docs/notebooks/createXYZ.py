f = open("BLBL.cart","r")
g = open("generate.xyz","w")

for i, line in enumerate (f):
   a = line.split()
   if i==0:
      nAt = a[3]
   elif (i==3):
      g.write(f"{a[1]} {a[2]} {a[3]}\n")
   elif (i==4):
      g.write(f"{a[1]} {a[2]} {a[3]}\n")
   elif (i==5):
      g.write(f"{a[1]} {a[2]} {a[3]}\n")
      g.write(f"{nAt}\n")
   elif (i>7):
      g.write(line)

f.close()
g.close()
    

