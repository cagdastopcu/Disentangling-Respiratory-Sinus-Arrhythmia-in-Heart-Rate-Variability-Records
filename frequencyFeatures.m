function [VLF, LF, HF, LFHFratio, nLF, nHF] = frequencyFeatures(f,PSD)      

      VLFInter = [0.0033 0.04];
      LFInter = [0.04 0.15];
      HFInter = [0.15 0.4];
      % find the indexes corresponding to the VLF, LF, and HF bands
      % cleaned with >VLFInter(1) zero Hz.
      iVLF= (f>=VLFInter(1)) & (f<=VLFInter(2));
      iLF = (f>=LFInter(1)) & (f<=LFInter(2));
      iHF = (f>=HFInter(1)) & (f<=HFInter(2));
      
      VLF=trapz(f(iVLF),PSD(iVLF));
      LF=trapz(f(iLF),PSD(iLF));
      HF=trapz(f(iHF),PSD(iHF));
      PSDtotal = trapz(f,PSD);
      LFHFratio = LF/HF; %sympathovagal balance 
      nLF = (LF/(PSDtotal-VLF))*100;
      nHF = (HF/(PSDtotal-VLF))*100;
end