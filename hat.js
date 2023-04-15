const r3 = 1.7320508075688772;
const hr3 = 0.8660254037844386;
const ident = [1,0,0,0,1,0];

let to_screen = [20, 0, 0, 0, -20, 0];
let lw_scale = 1;
let tiles;
let level;

let scale_centre;
let scale_start;
let scale_ts;

let reset_button;
let subst_button;
let translate_button;
let scale_button;
let draw_hats;
let draw_super;
let radio;
let cp_h1;
let cp_h;
let cp_t;
let cp_p;
let cp_f;
let dragging = false;
let uibox = true;
let box_height = 10;

let H1_hat;
let H_hat;
let T_hat;
let P_hat;
let F_hat;

let svg_serial = 0;
let h1_color = [0, 137, 212];
let h_color = [148, 205, 235];
let t_color = [251, 251, 251];
let p_color = [250, 250, 250];
let f_color = [191, 191, 191];

//-- hierachical coloring -- Begin
const PalRedCya = [[-1.0, 235, 20, 20], [0.0, 127, 127, 127], [1.0, 20, 255, 255]]; // default
const PalRedCyaInv=[[-1.0, 20, 255, 255], [0.0, 127, 127, 127], [1.0, 235, 20, 20]];
const PalGrey   = [[-1.0, 0, 0, 0], [1.0, 255, 255, 255]]; // geyscale
const PalGreyInv= [[-1.0, 255, 255, 255], [1.0, 0, 0, 0]]; // geyscale inverted
const PalGreRed = [ [-1.0, 215, 225, 238], [0.0, 134, 132, 138], [1.0, 153, 31, 23]]; // sequential Grey Red
const PalBluYel = [ [-1.0, 17, 95, 154], [0.0, 118, 198, 143], [1.0, 224, 241, 0]]; // classic sequential blue yellow 
const PalBlkPnk = [ [-1.0, 43, 43, 40], [0.0, 115, 110, 107], [1.0, 255, 128, 255]]; // sequential Black Pink
const PalBluRed = [ [-1.0, 25, 132, 197], [0.0, 226, 226, 226], [1.0, 194, 55, 40]]; // diverging Blue Red
const PalOrgPrp = [ [-1.0, 255, 180, 0], [0.0, 84, 82, 100], [1.0, 144, 128, 255]]; // diverging Orange Purple
const PaletteList = { PalDRedCya:PalRedCya, PalDRedCyaInv:PalRedCyaInv, PalSGrey:PalGrey, PalSGreyInv:PalGreyInv, PalSGreRed:PalGreRed, PalSBluYel:PalBluYel, PalSBlkPnk:PalBlkPnk, PalDBluRed:PalBluRed, PalDOrgPrp:PalOrgPrp };

const lv = -0.7;
const hv = 1.0;
const MotifH = { H1: hv/2, H: hv, T: lv, P: lv, F: lv };
const MotifT = { H1: lv, H: lv, T: hv, P: lv, F: lv };
const MotifP = { H1: lv, H: lv, T: lv, P: hv, F: lv };
const MotifF = { H1: lv, H: lv, T: lv, P: lv, F: hv };
const MotifPF = { H1: 0.0, H: 0.0, T: 0.0, P: 0.9, F: -0.9 };
const MotifHF = { H1: -0.4, H: 0.9, T: 0.0, P: 0.0, F: -0.9 };
const MotifDual = { H1: -0.7, H: 0.7, T: -0.0, P: -0.6, F: 0.6 };
const Motif0 = { H1: 0.0, H: 0.0, T: 0.0, P: 0.0, F: 0.0 }; // used for debug / check flat result
const MotifTest = { H1: 1.0, H: -1.0, T: -1.0, P: -0.5, F: 1.0 }; // test new motifs
const MotifList = { NoHierColor:null, MH:MotifH, MT:MotifT, MP:MotifP, MF:MotifF, MPF:MotifPF, MHF:MotifHF, MDual:MotifDual, M0:Motif0, Test:MotifTest };

let bHiCol = false; // enables Hierachical coloring
let Pal = PalRedCya; // selected Palette
let HatMotif = MotifDual; // selected motif for hat in cluster, range (-1,1)
let MetaMotif = MotifDual; // selected motif for metatiles, range (-1,1)

// adapts contrast to # levels and to Current Level, top metatile is neutral => #levels=level+1-1
function HiContrast(CurrLev = 0) {
  if (level == 1) return 1.0;
  const avgContrast = 1.0 / level; // linear normalization
  //const avgContrast = 1.0/sqrt(level); // for #levels high, normalization based on std_dev (2.0 bionomial)
  const prgContrast = -0.3; // range (-1: highlights top metatiles, 1: highlights hats), 0: neutral
  let LevContrast = (1 + prgContrast) - 2 * prgContrast * CurrLev / level;
  return avgContrast * LevContrast;
}
// maps val (range -1.0, 1.0) to color through palette Pal
function ColMap(val, Pal) {
  idx = max(min(val, 1.0), -1.0);
  i = 0;
  while (idx > Pal[i + 1][0]) i += 1;
  pc = (idx - Pal[i][0]) / (Pal[i + 1][0] - Pal[i][0]);
  HiFill = [
    Pal[i][1] + pc * (Pal[i + 1][1] - Pal[i][1]),
    Pal[i][2] + pc * (Pal[i + 1][2] - Pal[i][2]),
    Pal[i][3] + pc * (Pal[i + 1][3] - Pal[i][3]), ];
  return HiFill;
}
//-- hierachical coloring -- End

function color2rgb(c)
{
	return [c.levels[0], c.levels[1], c.levels[2]];
}
function set_h1_color()
{
	H1_hat.fill = color2rgb(cp_h1.color());
}
function set_h_color()
{
	H_hat.fill = color2rgb(cp_h.color());
}
function set_t_color()
{
	T_hat.fill = color2rgb(cp_t.color());
}
function set_p_color()
{
	P_hat.fill = color2rgb(cp_p.color());
}
function set_f_color()
{
	F_hat.fill = color2rgb(cp_f.color());
}

function pt( x, y )
{
	return { x : x, y : y };
}

function hexPt( x, y )
{
	return pt( x + 0.5*y, hr3*y );
}

// Affine matrix inverse
function inv( T ) {
	const det = T[0]*T[4] - T[1]*T[3];
	return [T[4]/det, -T[1]/det, (T[1]*T[5]-T[2]*T[4])/det,
		-T[3]/det, T[0]/det, (T[2]*T[3]-T[0]*T[5])/det];
};

// Affine matrix multiply
function mul( A, B )
{
	return [A[0]*B[0] + A[1]*B[3], 
		A[0]*B[1] + A[1]*B[4],
		A[0]*B[2] + A[1]*B[5] + A[2],

		A[3]*B[0] + A[4]*B[3], 
		A[3]*B[1] + A[4]*B[4],
		A[3]*B[2] + A[4]*B[5] + A[5]];
}

function padd( p, q )
{
	return { x : p.x + q.x, y : p.y + q.y };
}

function psub( p, q )
{
	return { x : p.x - q.x, y : p.y - q.y };
}

// Rotation matrix
function trot( ang )
{
	const c = cos( ang );
	const s = sin( ang );
    return [c, -s, 0, s, c, 0];
}

// Translation matrix
function ttrans( tx, ty )
{
    return [1, 0, tx, 0, 1, ty];
}

function rotAbout( p, ang )
{
	return mul( ttrans( p.x, p.y ), 
		mul( trot( ang ), ttrans( -p.x, -p.y ) ) );
}

// Matrix * point
function transPt( M, P )
{
	return pt(M[0]*P.x + M[1]*P.y + M[2], M[3]*P.x + M[4]*P.y + M[5]);
}

// Match unit interval to line segment p->q
function matchSeg( p, q )
{
	return [q.x-p.x, p.y-q.y, p.x,  q.y-p.y, q.x-p.x, p.y];
};

// Match line segment p1->q1 to line segment p2->q2
function matchTwo( p1, q1, p2, q2 )
{
	return mul( matchSeg( p2, q2 ), inv( matchSeg( p1, q1 ) ) );
};

function intersect( p1, q1, p2, q2 )
{
    const d = (q2.y - p2.y) * (q1.x - p1.x) - (q2.x - p2.x) * (q1.y - p1.y);
    const uA = ((q2.x - p2.x) * (p1.y - p2.y) - (q2.y - p2.y) * (p1.x - p2.x)) / d;
    const uB = ((q1.x - p1.x) * (p1.y - p2.y) - (q1.y - p1.y) * (p1.x - p2.x)) / d;

    return pt( p1.x + uA * (q1.x - p1.x), p1.y + uA * (q1.y - p1.y) );
}

function drawPolygon( shape, T, f, s, w )
{
  if( f != null ) {
	fill( ...f );
  } else {
    noFill();
  }
	if( s != null ) {
		stroke( ...s );
		strokeWeight( w * lw_scale );
  } else {
    noStroke();
  }
  beginShape();
	for( let p of shape ) {
		const tp = transPt( T, p );
		vertex( tp.x, tp.y );
  }
	endShape( CLOSE );
}

class Geom 
{
  constructor(pgon, fill, stroke, HiCol = 0) 
  {
    this.shape = pgon;
    this.fill = fill;
    this.stroke = stroke;
    this.width = 1.0;
    this.children = [];
    this.svg_id = null;
    this.HiCol = HiCol;
  }

  addChild( T, geom )
	{
		this.children.push( { T : T, geom : geom } );
  }

  evalChild( n, i )
  {
		return transPt( this.children[n].T, this.children[n].geom.shape[i] );
  }

  draw(S, level, HiCol = 0) {
    if (bHiCol) {
      if (level > 0) {
        for (let g of this.children) {
          g.geom.draw(mul(S, g.T), level - 1, this.HiCol * HiContrast(level) + HiCol);
        }
      } else {
        var Fill;
        if (this.fill == null) Fill = null;
        else Fill = ColMap(this.HiCol * HiContrast(level) + HiCol, Pal);
        drawPolygon(this.shape, S, Fill, this.stroke, this.width);
      }
    }
    else {
      if (level > 0) {
        for (let g of this.children) {
          g.geom.draw(mul(S, g.T), level - 1);
        }
      } else {
        drawPolygon( this.shape, S, this.fill, this.stroke, this.width );
      }
    }
  }

  recentre()
	{
    let cx = 0;
    let cy = 0;
	for( let p of this.shape ) {
      cx += p.x;
      cy += p.y;
    }
    cx /= this.shape.length;
    cy /= this.shape.length;
	const tr = pt( -cx, -cy );

	for( let idx = 0; idx < this.shape.length; ++idx ) {
		this.shape[idx] = padd( this.shape[idx], tr );
    }

	const M = ttrans( -cx, -cy );
	for( let ch of this.children ) {
		ch.T = mul( M, ch.T );
    }
  }

	resetSVG()
	{
	for( let ch of this.children ) {
      ch.geom.resetSVG();
    }
    this.svg_id = null;
    }

	buildSVGDefs( stream, sc )
	{
	if( this.svg_id != null ) {
      return;
    }

	this.svg_id = 't' + String(svg_serial).padStart( 5, '0' );
    ++svg_serial;

	for( let ch of this.children ) {
      const T = ch.T;
		ch.geom.buildSVGDefs( stream, sc * mag( T[0], T[1] ) );
    }

    // Stroked group
	stream.push( `  <g id="${this.svg_id}s">` );

	for( let ch of this.children ) {
      const T = ch.T;
      const gid = ch.geom.svg_id;
		stream.push( `    <use xlink:href="#${gid}s" transform="matrix(${T[0]} ${T[3]} ${T[1]} ${T[4]} ${T[2]} ${T[5]})"/>` );
    }

		if( (this.shape.length > 0) && (this.fill == null) ) {
			let verts = '';
			for( let p of this.shape ) {
				if( verts.length > 0 ) {
					verts = verts + ' ';
        }
				verts = verts + p.x + ',' + p.y;
      }
      let str = ' stroke="none"';
		if( this.stroke != null ) {
        let st = this.stroke;
				if( st.length == 1 ) {
          st = [st[0], st[0], st[0]];
        }

				 str = ` stroke="rgb(${st[0]},${st[1]},${st[2]})" stroke-width="${this.width*lw_scale/sc}"`;
      }

			stream.push( `    <polygon points="${verts}" fill="none"${str}/>` );
    }

		stream.push( '  </g>' );

    // Filled group
		stream.push( `  <g id="${this.svg_id}f">` );

		for( let ch of this.children ) {
      const T = ch.T;
      const gid = ch.geom.svg_id;
			stream.push( `    <use xlink:href="#${gid}f" transform="matrix(${T[0]} ${T[3]} ${T[1]} ${T[4]} ${T[2]} ${T[5]})"/>` );
    }

		if( (this.shape.length > 0) && (this.fill != null) ) {
			let verts = '';
			for( let p of this.shape ) {
				if( verts.length > 0 ) {
					verts = verts + ' ';
        }
				verts = verts + p.x + ',' + p.y;
      }
      let str = ' stroke="none"';
			if( this.stroke != null ) {
        let st = this.stroke;
				if( st.length == 1 ) {
          st = [st[0], st[0], st[0]];
        }

				 str = ` stroke="rgb(${st[0]},${st[1]},${st[2]})" stroke-width="${this.width*lw_scale/sc}"`;
      }

      let fil = ' fill="none"';
			if( this.fill != null ) {
        fil = ` fill="rgb(${this.fill[0]},${this.fill[1]},${this.fill[2]})"`;
      }

			stream.push( `    <polygon points="${verts}"${str}${fil}/>` );
    }

		stream.push( '  </g>' );
  }
}

const hat_outline = [
    hexPt(0, 0), hexPt(-1,-1), hexPt(0,-2), hexPt(2,-2),
    hexPt(2,-1), hexPt(4,-2), hexPt(5,-1), hexPt(4, 0),
    hexPt(3, 0), hexPt(2, 2), hexPt(0, 3), hexPt(0, 2),
    hexPt(-1, 2) ];

H1_hat = new Geom(hat_outline, h1_color, [0, 0, 0], HatMotif['H1']);
H_hat = new Geom(hat_outline, h_color, [0, 0, 0], HatMotif['H']);
T_hat = new Geom(hat_outline, t_color, [0, 0, 0], HatMotif['T']);
P_hat = new Geom(hat_outline, p_color, [0, 0, 0], HatMotif['P']);
F_hat = new Geom(hat_outline, f_color, [0, 0, 0], HatMotif['F']);

const H_init = (function () {
  const H_outline = [
		pt( 0, 0 ), pt( 4, 0 ), pt( 4.5, hr3 ),
		pt( 2.5, 5 * hr3 ), pt( 1.5, 5 * hr3 ), pt( -0.5, hr3 ) ];
  geom = new Geom(H_outline, null, [0, 0, 0], MetaMotif['H']);
  geom.width = 2;

  geom.addChild(
		matchTwo( 
			hat_outline[5], hat_outline[7], H_outline[5], H_outline[0] ),
		H_hat );
  geom.addChild(
		matchTwo( 
			hat_outline[9], hat_outline[11], H_outline[1], H_outline[2] ),
		H_hat );
  geom.addChild(
		matchTwo( 
			hat_outline[5], hat_outline[7], H_outline[3], H_outline[4] ),
		H_hat );
  geom.addChild(
		mul( ttrans( 2.5, hr3 ), 
    mul(
				[-0.5,-hr3,0,hr3,-0.5,0],
				[0.5,0,0,0,-0.5,0] ) ),
		H1_hat );

	return geom; }());

const T_init = (function () {
  const T_outline = [
		pt( 0, 0 ), pt( 3, 0 ), pt( 1.5, 3 * hr3 ) ];
  geom = new Geom(T_outline, null, [0, 0, 0], MetaMotif['T']);
  geom.width = 2;

	geom.addChild( 
		[0.5, 0, 0.5, 0, 0.5, hr3],
		T_hat );

	return geom; }());

const P_init = (function () {
	const P_outline = [
		pt( 0, 0 ), pt( 4, 0 ), 
		pt( 3, 2 * hr3 ), pt( -1, 2 * hr3 ) ];
  geom = new Geom(P_outline, null, [0, 0, 0], MetaMotif['P']);
  geom.width = 2;

  geom.addChild(
		[0.5, 0, 1.5, 0, 0.5, hr3],
		P_hat );
	geom.addChild( 
		mul( ttrans( 0, 2 * hr3 ), 
			mul( [0.5, hr3, 0, -hr3, 0.5, 0],
				 [0.5, 0.0, 0.0, 0.0, 0.5, 0.0] ) ),
		P_hat );

	return geom; }());

const F_init = (function () {
  const F_outline = [
		pt( 0, 0 ), pt( 3, 0 ), 
		pt( 3.5, hr3 ), pt( 3, 2 * hr3 ), pt( -1, 2 * hr3 ) ];
  geom = new Geom(F_outline, null, [0, 0, 0], MetaMotif['F']);
  geom.width = 2;

  geom.addChild(
		[0.5, 0, 1.5, 0, 0.5, hr3],
		F_hat );
	geom.addChild( 
		mul( ttrans( 0, 2 * hr3 ), 
			mul( [0.5, hr3, 0, -hr3, 0.5, 0],
				 [0.5, 0.0, 0.0, 0.0, 0.5, 0.0] ) ),
		F_hat );

	return geom; }());

function constructPatch( H, T, P, F )
{
  const rules = [
		['H'],
		[0, 0, 'P', 2],
		[1, 0, 'H', 2],
		[2, 0, 'P', 2],
		[3, 0, 'H', 2],
		[4, 4, 'P', 2],
		[0, 4, 'F', 3],
		[2, 4, 'F', 3],
		[4, 1, 3, 2, 'F', 0],
		[8, 3, 'H', 0],
		[9, 2, 'P', 0],
		[10, 2, 'H', 0],
		[11, 4, 'P', 2],
		[12, 0, 'H', 2],
		[13, 0, 'F', 3],
		[14, 2, 'F', 1],
		[15, 3, 'H', 4],
		[8, 2, 'F', 1], 
		[17, 3, 'H', 0],
		[18, 2, 'P', 0],
		[19, 2, 'H', 2],
		[20, 4, 'F', 3],
		[20, 0, 'P', 2],
		[22, 0, 'H', 2],
		[23, 4, 'F', 3],
		[23, 0, 'F', 3],
		[16, 0, 'P', 2],
		[9, 4, 0, 2, 'T', 2],
		[4, 0, 'F', 3] 
  ];

	ret = new Geom( [], null, null );
  ret.width = H.width;
	shapes = { 'H' : H, 'T' : T, 'P' : P, 'F' : F };

	for( let r of rules ) {
		if( r.length == 1 ) {
			ret.addChild( ident, shapes[r[0]] );
		} else if( r.length == 4 ) {
      const poly = ret.children[r[0]].geom.shape;
      const T = ret.children[r[0]].T;
			const P = transPt( T, poly[(r[1]+1)%poly.length] );
			const Q = transPt( T, poly[r[1]] );
      const nshp = shapes[r[2]];
      const npoly = nshp.shape;

      ret.addChild(
				matchTwo( npoly[r[3]], npoly[(r[3]+1)%npoly.length], P, Q ),
				nshp );
    } else {
      const chP = ret.children[r[0]];
      const chQ = ret.children[r[2]];

	  const P = transPt( chQ.T, chQ.geom.shape[r[3]] );
	  const Q = transPt( chP.T, chP.geom.shape[r[1]] );
      const nshp = shapes[r[4]];
      const npoly = nshp.shape;

      ret.addChild(
				matchTwo( npoly[r[5]], npoly[(r[5]+1)%npoly.length], P, Q ),
				nshp );
    }
  }

  return ret;
}

function constructMetatiles( patch )
{
	const bps1 = patch.evalChild( 8, 2 );
	const bps2 = patch.evalChild( 21, 2 );
	const rbps = transPt( rotAbout( bps1, -2.0*PI/3.0 ), bps2 );

	const p72 = patch.evalChild( 7, 2 );
	const p252 = patch.evalChild( 25, 2 );

	const llc = intersect( bps1, rbps,
	patch.evalChild( 6, 2 ), p72 );
	let w = psub( patch.evalChild( 6, 2 ), llc );

    const new_H_outline = [llc, bps1];
	w = transPt( trot( -PI/3 ), w );
	new_H_outline.push( padd( new_H_outline[1], w ) );
	new_H_outline.push( patch.evalChild( 14, 2 ) );
	w = transPt( trot( -PI/3 ), w );
	new_H_outline.push( psub( new_H_outline[3], w ) );
	new_H_outline.push( patch.evalChild( 6, 2 ) );

    const new_H = new Geom(new_H_outline, null, [0, 0, 0], MetaMotif['H']);
    new_H.width = patch.width * 2;
	for( let ch of [0, 9, 16, 27, 26, 6, 1, 8, 10, 15] ) {
		new_H.addChild( patch.children[ch].T, patch.children[ch].geom );
    }

	const new_P_outline = [ p72, padd( p72, psub( bps1, llc ) ), bps1, llc ];
    const new_P = new Geom(new_P_outline, null, [0, 0, 0], MetaMotif['P']);
    new_P.width = patch.width * 2;
	for( let ch of [7,2,3,4,28] ) {
		new_P.addChild( patch.children[ch].T, patch.children[ch].geom );
    }

    const new_F_outline = [
		bps2, patch.evalChild( 24, 2 ), patch.evalChild( 25, 0 ),
		p252, padd( p252, psub( llc, bps1 ) ) ];
    const new_F = new Geom(new_F_outline, null, [0, 0, 0], MetaMotif['F']);
    new_F.width = patch.width * 2;
	for( let ch of [21,20,22,23,24,25] ) {
		new_F.addChild( patch.children[ch].T, patch.children[ch].geom );
    }

    const AAA = new_H_outline[2];
	const BBB = padd( new_H_outline[1], 
	psub( new_H_outline[4], new_H_outline[5] ) );
	const CCC = transPt( rotAbout( BBB, -PI/3 ), AAA );
	const new_T_outline = [BBB,CCC,AAA];
    const new_T = new Geom(new_T_outline, null, [0, 0, 0], MetaMotif['T']);
    new_T.width = patch.width * 2;
	new_T.addChild( patch.children[11].T, patch.children[11].geom );

    new_H.recentre();
    new_P.recentre();
    new_F.recentre();
    new_T.recentre();

	return [new_H, new_T, new_P, new_F]
}

function isButtonActive( but )
{
  return but.elt.style.border.length > 0;
}

function setButtonActive( but, b )
{
	but.elt.style.border = (b ? "3px solid black" : "");
}

function setup() {
	createCanvas( windowWidth, windowHeight );

    tiles = [H_init, T_init, P_init, F_init];
    level = 1;

	reset_button = createButton( "Reset" );
	reset_button.position( 10, box_height );
	reset_button.size( 125, 25 );
	reset_button.mousePressed( function() {
    tiles = [H_init, T_init, P_init, F_init];
    level = 1;
		radio.selected( 'H' );
    to_screen = [20, 0, 0, 0, -20, 0];
    lw_scale = 1;
		setButtonActive( draw_hats, true );
		setButtonActive( draw_super, true );
    loop();
	} );
	box_height += 30;

	subst_button = createButton( "Build Supertiles" );
	subst_button.position( 10, box_height );
	subst_button.size( 125, 25 );
	subst_button.mousePressed( function() {
		const patch = constructPatch( ...tiles );
		tiles = constructMetatiles( patch );
    ++level;
    loop();
	} );
	box_height += 40;

    radio = createRadio();
	radio.mousePressed( function() { loop() } );
	radio.position( 10, box_height );
	for( let s of ['H', 'T', 'P', 'F'] ) {
		let o = radio.option( s );
    o.onclick = loop;
    }
	radio.selected( 'H' );
	box_height += 40;

	cp_h1 = createColorPicker( color( ...h1_color ) );
	cp_h1.mousePressed( function() { loop() } );
	cp_h1.input(set_h1_color);
	cp_h1.position( 10, box_height );
	box_height += 40;

	cp_h = createColorPicker( color( ...h_color) );
	cp_h.mousePressed( function() { loop() } );
	cp_h.input(set_h_color);
	cp_h.position( 10, box_height );
	box_height += 40;

	cp_t = createColorPicker( color( ...t_color ) );
	cp_t.mousePressed( function() { loop() } );
	cp_t.input(set_t_color);
	cp_t.position( 10, box_height );
	box_height += 40;

	cp_p = createColorPicker( color( ...p_color) );
	cp_p.mousePressed( function() { loop() } );
	cp_p.input(set_p_color);
	cp_p.position( 10, box_height );
	box_height += 40;

	cp_f = createColorPicker( color( ...f_color ) );
	cp_f.mousePressed( function() { loop() } );
	cp_f.input(set_f_color);
	cp_f.position( 10, box_height );
	box_height += 40;

	translate_button = createButton( "Translate" );
	setButtonActive( translate_button, true );
	translate_button.position( 10, box_height );
	translate_button.size( 125, 25 );
	translate_button.mousePressed( function() {
		setButtonActive( translate_button, true );
		setButtonActive( scale_button, false );
    loop();
	} );
	box_height += 30;

	scale_button = createButton( "Scale" );
	scale_button.position( 10, box_height );
	scale_button.size( 125, 25 );
	scale_button.mousePressed( function() {
		setButtonActive( translate_button, false );
		setButtonActive( scale_button, true );
    loop();
	} );
	box_height += 30;

	draw_hats = createButton( "Draw Hats" );
	setButtonActive( draw_hats, true );
	draw_hats.mousePressed( function() {
		setButtonActive( draw_hats, !isButtonActive( draw_hats ) );
    loop();
	} );
	draw_hats.position( 10, box_height );
	draw_hats.size( 125, 25 );
	box_height += 30;

	draw_super = createButton( "Draw Supertiles" );
	setButtonActive( draw_super, true );
	draw_super.mousePressed( function() {
		setButtonActive( draw_super, !isButtonActive( draw_super ) );
    loop();
	} );
	draw_super.position( 10, box_height );
	draw_super.size( 125, 25 );
	box_height += 40;

	let save_button = createButton( "Save PNG" );
	save_button.position( 10, box_height );
	save_button.size( 125, 25 );
	save_button.mousePressed( function () {
    uibox = false;
    draw();
    if (bHiCol) save("Hat_L" + str(level) + "_" + radio.value() + "_" + selMotif.value() + "_" + selPalette.value() + ".png");
    else save("Hat_L" + str(level) + "_" + radio.value() + "_byCluster" + ".png");
    uibox = true;
    draw();
	} );
	box_height += 30;

	let svg_button = createButton( "Save SVG" );
	svg_button.position( 10, box_height );
	svg_button.size( 125, 25 );
	svg_button.mousePressed( function () {
    svg_serial = 0;
		for( let t of tiles ) {
      t.resetSVG();
    }

    const stream = [];
		stream.push( `<svg viewBox="0 0 ${width} ${height}" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">` );
		stream.push( '<defs>' );
		for( let t of tiles ) {
			t.buildSVGDefs( stream, mag( to_screen[0], to_screen[1] ) );
		}
		stream.push( '</defs>' );

		const idx = {'H':0, 'T':1, 'P':2, 'F':3}[radio.value()];

		if( isButtonActive( draw_hats ) ) {
			stream.push( `<use xlink:href="#${tiles[idx].svg_id}f" transform="matrix(${to_screen[0]} ${to_screen[3]} ${to_screen[1]} ${to_screen[4]} ${to_screen[2]+width/2} ${to_screen[5]+height/2})"/>` )
		}
		if( isButtonActive( draw_super ) ) {
			stream.push( `<use xlink:href="#${tiles[idx].svg_id}s" transform="matrix(${to_screen[0]} ${to_screen[3]} ${to_screen[1]} ${to_screen[4]} ${to_screen[2]+width/2} ${to_screen[5]+height/2})"/>` )
		}
		stream.push( '</svg>' );

		saveStrings( stream, 'output', 'svg' );
	} );
	box_height += 30;
	box_height -= 5; // remove half the padding
  
  //-- Hierchical color options --
  let selMotif = createSelect();
  selMotif.position(10, 550);
  selMotif.size(125, 20);
  selMotif.selected(Object.keys(MotifList)[0]);
  for (let a of Object.keys(MotifList)) selMotif.option(a);
  selMotif.changed(function mySelectEvent() {
    let item = selMotif.value();
    if (item == Object.keys(MotifList)[0]) bHiCol = false;
    else {
      bHiCol = true;
      MetaMotif = MotifList[item];
      HatMotif = MotifList[item]; // for future option : manage a different motif for Lev 0
	  H1_hat.HiCol = MetaMotif['H1'];
      H_hat.HiCol = MetaMotif['H'];
      T_hat.HiCol = MetaMotif['T'];
      P_hat.HiCol = MetaMotif['P'];
      F_hat.HiCol = MetaMotif['F'];
      tiles = [H_init, T_init, P_init, F_init];
      tiles[0].HiCol = MetaMotif['H'];
      tiles[1].HiCol = MetaMotif['T'];
      tiles[2].HiCol = MetaMotif['P'];
      tiles[3].HiCol = MetaMotif['F'];
      levelnow = level;
      level = 1;
      for (let lev = levelnow; lev > 1; --lev) {
        const patch = constructPatch(...tiles);
        tiles = constructMetatiles(patch);
        ++level;
      }
    }
    loop();
  });

  let selPalette = createSelect();
  selPalette.position(10, 580);
  selPalette.size(125, 20);
  selPalette.selected(Object.keys(PaletteList)[0]);
  for (let a of Object.keys(PaletteList)) selPalette.option(a);
  selPalette.changed(function mySelectEvent() {
    let item = selPalette.value();
    Pal = PaletteList[item];
    loop();
  });
}

function draw()
{
    background( 255 );

    push();
	translate( width/2, height/2 );
	const idx = {'H':0, 'T':1, 'P':2, 'F':3}[radio.value()];
    tiles[idx].HiCol = 0; // top metatile is neutral in hierchical coloring

	if( isButtonActive( draw_hats ) ) {
		tiles[idx].draw( to_screen, level );
    }

	if( isButtonActive( draw_super ) ) {
		for( let lev = level - 1; lev >= 0; --lev ) {
			tiles[idx].draw( to_screen, lev );
        }
    }
    pop();

	if( uibox ) {
		stroke( 0 );
		strokeWeight( 0.5 );
		fill( 255, 220 );
		rect( 5, 5, 135, box_height);
    }
   noLoop();
}

function windowResized() 
{
	resizeCanvas( windowWidth, windowHeight );
}

function mousePressed()
{
  dragging = true;
	if( isButtonActive( scale_button ) ) {
		scale_centre = transPt( inv( to_screen ), pt( width/2, height/2 ) );
		scale_start = pt( mouseX, mouseY );
    scale_ts = [...to_screen];
  }
  loop();
}

function mouseDragged()
{
	if( dragging ) {
		if( isButtonActive( translate_button ) ) {
			to_screen = mul( ttrans( mouseX - pmouseX, mouseY - pmouseY ), 
				to_screen );
		} else if( isButtonActive( scale_button ) ) {
			let sc = dist( mouseX, mouseY, width/2, height/2 ) / 
				dist( scale_start.x, scale_start.y, width/2, height/2 );
      to_screen = mul(
				mul( ttrans( scale_centre.x, scale_centre.y ),
					mul( [sc, 0, 0, 0, sc, 0],
						ttrans( -scale_centre.x, -scale_centre.y ) ) ),
				scale_ts );
			lw_scale = mag( to_screen[0], to_screen[1] ) / 20.0;
    }
    loop();
    return false;
  }
}

function mouseReleased()
{
  dragging = false;
  //loop(); // original sw commented, this line seems not necessary? it produces additional refreshes
}
