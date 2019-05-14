/**
 * git clone https://github.com/hou-dao/deom.git
 * ---
 * Written by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include "deom.hpp"

void deom::oprAct (cx_cube& d_ddos, const cx_mat& sdip, const cx_cube& ddos, const char lrc) 
{
    for (int iado=0; iado<nddo; ++iado) 
    {
        const cx_mat& ado = ddos.slice(iado);
        if (lrc == 'l') 
        {
            d_ddos.slice(iado) += sdip*ado;
        } 
        else if (lrc == 'r') 
        {
            d_ddos.slice(iado) += ado*sdip;
        } 
        else if (lrc == 'c') 
        {
            d_ddos.slice(iado) += sdip*ado-ado*sdip;
        }
    }
}

void deom::iniHei (cx_cube& ddos, const cx_mat& sdip)
{
    ddos.slice(0) = sdip;
    nddo = 1;
}