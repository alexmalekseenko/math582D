% returns Legendre polynomial of degree n
function p = legendre_p ( n )

    % for now just choose from known polynomials
    switch n
        case 0
            p = [1] ;
        case 1
            p = [1 0] ;
        case 2
            p = [3 0 -1]/2 ;
        case 3
            p = [5 0 -3 0]/2 ;
        case 4
            p = [35 0 -30 0 3]/8 ;
        case 5
            p = [63 0 -70 0 15 0]/8 ;
        case 6
            p = [231 0 -315 0 105 0 -5]/16 ;
		case 7
			p = [429 0 -693 0 315 0 -35 0]/16 ;
		case 8
			p = [6435 0 -12012 0 6930 0 -1260 0 35]/128 ;
		otherwise
			error('Unsupported Legendre Polynomial degree') ;
    end


