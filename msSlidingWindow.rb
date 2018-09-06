#!/usr/bin/env ruby
#
#

load "~/rubyStuff/Hudson.rb"


#bring in file
m = MSRun.new(ARGV[0])
m.samples.each{|aSample|
  w = 0.0
  i = 0
  while(w<1.0)
    count = Array.new(aSample.seqMat.sampleSize - 1,0)

    pos = aSample.positions.select{|x| x.to_f > w && x.to_f < w+0.1}
    pos.each{ | aPos |
      ind= aSample.positions.index(aPos)
      1.upto(aSample.seqMat.sampleSize-1){ | j |
#        print aSample.seqMat.matrix[1][1,1]
        if (aSample.seqMat.matrix[0][ind,1] != aSample.seqMat.matrix[j][ind,1])
          count[j-1] = count[j-1].to_i + 1
        end
      }
    }
    w+=0.05
    count.each{ |anAllele| print anAllele,"\t" }
    print "\n"
  end
}