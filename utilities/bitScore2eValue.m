function eValue = bitScore2eValue(thres, qSize, DbSize)
  eValue = qSize*DbSize*2^(-thres);
end

