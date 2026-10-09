f1 = open("BLBL.cart", "r")
f2 = open("BLBL.mol", "r")
g = open("BLBLNew.mol", "w")

for i in range(8):
   f1.readline()

for line in f2:
   a = line.split()
   if (len(a) == 10):
      b = f1.readline().split()
      g.write("{} {} {} {} {} {} {} {} {} {}\n".format(a[0], a[1], a[2], a[3], b[1], b[2], b[3], a[7], a[8], a[9]))
   else:
      g.write(line)
