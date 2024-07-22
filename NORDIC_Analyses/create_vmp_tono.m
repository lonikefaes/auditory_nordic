function  create_vmp_tono(filenamepath,NameFile,VTCInfo,values,index,mapname)

dircolor = 'C:\Program Files\BrainVoyager21\MapLUTs\jet_posneg2.olt';

vmp = xff('vmp');
vmp.XStart = VTCInfo.XStart;
vmp.XEnd = VTCInfo.XEnd;
vmp.YStart = VTCInfo.YStart;
vmp.YEnd = VTCInfo.YEnd;
vmp.ZStart = VTCInfo.ZStart;
vmp.ZEnd = VTCInfo.ZEnd;
vmp.Resolution = VTCInfo.Resolution;

    for it = 1:size(values,2)
        if it>1
            vmp.Map(it) = vmp.Map(it-1);
        end
        vmp.Map(it).VMPData              = zeros(VTCInfo.size);
        vmp.Map(it).VMPData(index) = values(:,it);
        vmp.Map(it).LUTName = dircolor;   

        if(exist('mapname'))
        vmp.Map(it).Name                 = mapname{it};     
        end
    end
    
vmp.NrOfMaps = size(values,2);
vmp.SaveAs([filenamepath,NameFile]);
vmp.clearobject;
clear vmp
