function bdry = generateRotGaussBoundary( g )

K=g.numT;
bdry = false( K, 3 );

for iT = 1:K
    for iE=1:3
        %If not interior edge
       if (g.idE0T(iT, iE) ~= 0) 
           edgeNr = g.E0T(iT, iE);
           
           baryX = g.baryE( edgeNr, 1 );
           baryY = g.baryE( edgeNr, 2 );
           
           %Left
           if ( baryX == -0.5 && baryY > 0.0)
               bdry( iT, iE ) = 1;
           end
           %Right
           if ( baryX ==  0.5 && baryY < 0.0)
               bdry( iT, iE ) = 1;
           end
           %Bottom
           if ( baryY == -0.5 && baryX < 0.0)
               bdry( iT, iE ) = 1;
           end
           %Top
           if ( baryY ==  0.5 && baryX > 0.0)
               bdry( iT, iE ) = 1;
           end
           
       end
    end
    
end

end