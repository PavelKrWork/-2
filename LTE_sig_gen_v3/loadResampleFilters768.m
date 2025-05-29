function IRFFilterSubBand = loadResampleFilters768()
              
%         load( 'FilterCoef768\BBSubBand2Filter1FloatPointI1D1.mat' );
%         
%         IRFFilterSubBand{ 1, 1} = FloatPointResampleFilter;

        IRFFilterSubBand{ 1, 1} = 1.0;
        
        load( 'LTE_sig_gen_v3\FilterCoef768\IRFSubBand1Filter1FloatPointI4D3.mat' );
        
        IRFFilterSubBand{ 1, 2} = FloatPointResampleFilter;
        
        load( 'LTE_sig_gen_v3\FilterCoef768\IRFSubBand1Filter2FloatPointI2D1.mat' );
        
        IRFFilterSubBand{ 1, 3} = FloatPointResampleFilter;

        load( 'LTE_sig_gen_v3\FilterCoef768\IRFSubBand1Filter3FloatPointI2D1.mat' );
        
        IRFFilterSubBand{ 1, 4} = FloatPointResampleFilter;
               
end