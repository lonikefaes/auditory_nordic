function VTCinfo = extractVTCinfo(vtc)
    VTCinfo.XStart       = vtc.XStart; VTCinfo.XEnd      = vtc.XEnd;  
    VTCinfo.YStart       = vtc.YStart; VTCinfo.YEnd = vtc.YEnd; 
    VTCinfo.ZStart       = vtc.ZStart; VTCinfo.ZEnd   = vtc.ZEnd;
    VTCinfo.Resolution   = vtc.Resolution; temp = size(vtc.VTCData);VTCinfo.size =  temp(2:end);
    VTCinfo.TR           = vtc.TR;
    VTCinfo.NrOfVolumes  = size(vtc.VTCData,1);
    