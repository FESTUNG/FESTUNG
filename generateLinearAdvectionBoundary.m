function bdry = generateLinearAdvectionBoundary( g )

K=g.numT;
bdry = false( K, 3 );

for iT = 1:K
    for iE=1:3
        %If not interior edge
       if (g.idE0T(iT, iE) ~= 0) 
           edgeNr = g.E0T(iT, iE);
           
           baryX = g.baryE( edgeNr, 1 );
           baryY = g.baryE( edgeNr, 2 );
           
           %Right
           if ( baryX > 0.4999999 )
               bdry( iT, iE ) = 1;
           end
           %Top
           if ( baryY > 0.4999999 )
               bdry( iT, iE ) = 1;
           end
 
       end
    end
    
end

end