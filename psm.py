import cv2
import numpy as np
import matplotlib.pyplot as plt
import os



class PSM():
   
    #Could be optimised by creating submethods for most parts between comments
    
    def psm(self,dossier): 
        entries = os.listdir(dossier) 
        entries_filtered = []
        for p in entries:
            entries_filtered.append(int(p.replace('.png', '')))
            
        entries_filtered.sort()

        stop = np.amax(entries_filtered)
        start = np.amin(entries_filtered)
        step = entries_filtered[1]-entries_filtered[0]
        
        list_images = []
        for k in entries_filtered:
            if k > 0:
                list_images.append('+'+str(k)+'.png')
            else:
                list_images.append(str(k)+'.png')

        l = np.size(list_images) 
        f, axarr = plt.subplots(7,l) # Create the figure
        

        taille = [] #Init some variables for future use.
        vmin = []
        vmax= []
        j = 0
        norm = 0
        pui = 0
        lp = []
        
        for s in list_images: # Find max value for normalisation
            ima = dossier+'/'+s
            immax = cv2.imread(ima, cv2.IMREAD_GRAYSCALE)
            if norm < np.amax(immax):
                norm = np.amax(immax)
            r = np.amax(immax)
        ta = []
        for k in list_images: # First instance with guessed threshold
            ima = dossier+'/'+k
            immax = cv2.imread(ima, cv2.IMREAD_GRAYSCALE)
            r = np.amax(immax)
            rett, thresht = cv2.threshold(immax, np.amax(immax)/2, 255, 0)
            Mt = cv2.moments(thresht)  
            if Mt["m00"] != 0:
                cX = int(Mt["m10"] / Mt["m00"])
                cY = int(Mt["m01"] / Mt["m00"])
            else:
                cX, cY = 0, 0   
            y=cY
            x=cX
            
            immax[immax < 5] = 0
            somme = np.sum(immax**2)
            for wt in range(20,100,2): # Estimate the cropped image size
                test = immax[int(y-wt/2):int(y+wt/2)+1, int(x-wt/2):int(x+wt/2)+1]

                taille.append(np.sum(test**2))
            w = np.where(taille > 0.85*somme)[0][0]*2+20
            ta.append(w)
        
        w = int(np.mean(ta))
    
        
            
        
        
            

        for s in list_images:
            ima = dossier+'/'+s
            immax = cv2.imread(ima, cv2.IMREAD_GRAYSCALE)
            r = np.amax(immax)
            rett, thresht = cv2.threshold(immax, np.amax(immax)/2, 255, 0)
            Mt = cv2.moments(thresht)  
            if Mt["m00"] != 0:
                cX = int(Mt["m10"] / Mt["m00"])
                cY = int(Mt["m01"] / Mt["m00"])
            else:
                cX, cY = 0, 0   
            y=cY
            x=cX
                
                
            
            test = immax[int(y-w/2):int(y+w/2)+1, int(x-w/2):int(x+w/2)+1] #Base power present in the images for further analysis
            test = test/norm
            if pui < np.sum(test**2):
                pui = np.sum(test**2)
            lp.append(pui)
            pui = 0
        
        y1 = []
        y2 = []
        y3 = []
        y4 = []
        y5 = []
        y6 = []
        y7 = []

        
        abb = []
        res = []
        

        for i in list_images:
            for c in range(1, r): # Test for best threshold
                pic = dossier+'/'+i
                im = cv2.imread(pic, 0)
                threshold = c
                ret, thresh = cv2.threshold(im, threshold, 255, 0)
                M = cv2.moments(thresh)

        
                
                if M["m00"] != 0:
                    cX = int(M["m10"] / M["m00"])
                    cY = int(M["m01"] / M["m00"])
                else:
                    cX, cY = 0, 0   
                y=cY
                x=cX
                
                
            
                raw= im[int(y-w/2):int(y+w/2)+1, int(x-w/2):int(x+w/2)+1]
                

                raw = raw/norm
                rot1 = cv2.rotate(raw, cv2.cv2.ROTATE_90_CLOCKWISE)
                rot2 = cv2.rotate(rot1, cv2.cv2.ROTATE_90_CLOCKWISE)
                rot3 = cv2.rotate(rot2, cv2.cv2.ROTATE_90_CLOCKWISE)

                sym = (raw+rot1+rot2+rot3)/4
            
                asym = raw-sym


                Fee = (asym + np.fliplr(asym) + np.flipud(asym)+ np.fliplr(np.flipud(asym)))/4
                Foo = (asym - np.fliplr(asym) - np.flipud(asym) + np.fliplr(np.flipud(asym)))/4
                Feo = (asym + np.fliplr(asym) - np.flipud(asym) - np.fliplr(np.flipud(asym)))/4
                Foe = (asym - np.fliplr(asym) + np.flipud(asym) - np.fliplr(np.flipud(asym)))/4


                if np.sum(raw**2,dtype=np.float64) <  lp[j]/1.25: # Verification that most of the power is present in the frame to confirm a good centering
                    abb.append(100)
                elif (np.sum(Feo**2,dtype=np.float64))+(np.sum(Foe**2,dtype=np.float64)) < 0.2:
                    abb.append(101)
                else:
                    abb.append((np.sum(Feo**2,dtype=np.float64))+(np.sum(Foe**2,dtype=np.float64))) 

            res.append(np.argmin(abb)+1)
            abb = []

        
        ceny = 100
        cenx = 100
        centery =[]
        centerx =[]
        for i in list_images: # Test for a few point around the current centroid
            pic = dossier+'/'+i
            im = cv2.imread(pic, 0)
            threshold = res[j]
            ret, thresh = cv2.threshold(im, threshold, 255, 0)
            M = cv2.moments(thresh)

    
            
            if M["m00"] != 0:
                cX = int(M["m10"] / M["m00"])
                cY = int(M["m01"] / M["m00"])
            else:
                cX, cY = 0, 0   
            y=cY
            x=cX
            
            for q in range(-int(w/5),int(w/5)): # Find best centroid in y
                y = y+q
                raw= im[int(y-w/2):int(y+w/2)+1, int(x-w/2):int(x+w/2)+1]
                

                raw = raw/norm
                rot1 = cv2.rotate(raw, cv2.cv2.ROTATE_90_CLOCKWISE)
                rot2 = cv2.rotate(rot1, cv2.cv2.ROTATE_90_CLOCKWISE)
                rot3 = cv2.rotate(rot2, cv2.cv2.ROTATE_90_CLOCKWISE)

                sym = (raw+rot1+rot2+rot3)/4
                
                asym = raw-sym


                Fee = (asym + np.fliplr(asym) + np.flipud(asym)+ np.fliplr(np.flipud(asym)))/4
                Foo = (asym - np.fliplr(asym) - np.flipud(asym) + np.fliplr(np.flipud(asym)))/4
                Feo = (asym + np.fliplr(asym) - np.flipud(asym) - np.fliplr(np.flipud(asym)))/4
                Foe = (asym - np.fliplr(asym) + np.flipud(asym) - np.fliplr(np.flipud(asym)))/4

                
                    
                if ceny > (np.sum(Feo**2,dtype=np.float64)):
                    ceny = (np.sum(Feo**2,dtype=np.float64))
                    resy = y
                    y = cY
                else:
                    y = cY
            ceny = 100
            centery.append(resy)
            for h in range(-int(w/5),int(w/5)): # Find best centroid in x
                x = x+h
                raw= im[int(y-w/2):int(y+w/2)+1, int(x-w/2):int(x+w/2)+1]
                

                raw = raw/norm
                rot1 = cv2.rotate(raw, cv2.cv2.ROTATE_90_CLOCKWISE)
                rot2 = cv2.rotate(rot1, cv2.cv2.ROTATE_90_CLOCKWISE)
                rot3 = cv2.rotate(rot2, cv2.cv2.ROTATE_90_CLOCKWISE)

                sym = (raw+rot1+rot2+rot3)/4
                
                asym = raw-sym


                Fee = (asym + np.fliplr(asym) + np.flipud(asym)+ np.fliplr(np.flipud(asym)))/4
                Foo = (asym - np.fliplr(asym) - np.flipud(asym) + np.fliplr(np.flipud(asym)))/4
                Feo = (asym + np.fliplr(asym) - np.flipud(asym) - np.fliplr(np.flipud(asym)))/4
                Foe = (asym - np.fliplr(asym) + np.flipud(asym) - np.fliplr(np.flipud(asym)))/4

                    
                if cenx > (np.sum(Foe**2,dtype=np.float64)):
                    cenx = (np.sum(Foe**2,dtype=np.float64))
                    resx = x
                    x = cX
                else:
                    x = cX
            cenx = 100 
            centerx.append(resx)
        
        
        for i in list_images: # Fills the figures
            pic = dossier+'/'+i
            im = cv2.imread(pic, 0)
            
            y = centery[j]
            x = centerx[j]
            
            
            raw= im[int(y-w/2):int(y+w/2)+1, int(x-w/2):int(x+w/2)+1]
            

            raw = raw/norm
            rot1 = cv2.rotate(raw, cv2.cv2.ROTATE_90_CLOCKWISE)
            rot2 = cv2.rotate(rot1, cv2.cv2.ROTATE_90_CLOCKWISE)
            rot3 = cv2.rotate(rot2, cv2.cv2.ROTATE_90_CLOCKWISE)

            sym = (raw+rot1+rot2+rot3)/4
        
            asym = raw-sym


            Fee = (asym + np.fliplr(asym) + np.flipud(asym)+ np.fliplr(np.flipud(asym)))/4
            Foo = (asym - np.fliplr(asym) - np.flipud(asym) + np.fliplr(np.flipud(asym)))/4
            Feo = (asym + np.fliplr(asym) - np.flipud(asym) - np.fliplr(np.flipud(asym)))/4
            Foe = (asym - np.fliplr(asym) + np.flipud(asym) - np.fliplr(np.flipud(asym)))/4

            
            
            y1.append(np.sum(raw**2,dtype=np.float64))
            y2.append(np.sum(sym**2,dtype=np.float64))
            y3.append(np.sum(asym**2,dtype=np.float64))
            y4.append(np.sum(Fee**2,dtype=np.float64))
            y5.append(np.sum(Foo**2,dtype=np.float64))
            y6.append(np.sum(Feo**2,dtype=np.float64))
            y7.append(np.sum(Foe**2,dtype=np.float64))

            vmin.append(raw.min())
            vmax.append(raw.max())
            vmin.append(sym.min())
            vmax.append(sym.max())
            vmin.append(asym.min())
            vmax.append(asym.max())
            vmin.append(Fee.min())
            vmax.append(Fee.max())
            vmin.append(Foo.min())
            vmax.append(Foo.max())
            vmin.append(Feo.min())
            vmax.append(Feo.max())
            vmin.append(Foe.min())
            vmax.append(Foe.max())
            

            nom = ['raw', 'sym', 'asym', 'ee', 'oo', 'eo', 'oe']

            
            axarr[0][j].imshow(raw, cmap='jet',vmin=0,vmax=1)
            axarr[0][j].axis('off')
            axarr[0][j].text(0.5, 1.1, i.replace('.png','$\;\mu m$'),horizontalalignment='center', transform=axarr[0][j].transAxes)

            axarr[1][j].imshow(sym, cmap='jet',vmin=0,vmax=1)
            axarr[1][j].axis('off')

            axarr[2][j].imshow(asym, cmap='jet', vmin = vmin[2], vmax = vmax[2])
            axarr[2][j].axis('off')

            axarr[3][j].imshow(Fee, cmap='jet', vmin = vmin[2], vmax = vmax[2])
            axarr[3][j].axis('off')

            axarr[4][j].imshow(Foo, cmap='jet', vmin = vmin[2], vmax = vmax[2])
            axarr[4][j].axis('off')

            axarr[5][j].imshow(Feo, cmap='jet', vmin = vmin[2], vmax = vmax[2])
            axarr[5][j].axis('off')
            

            axarr[6][j].imshow(Foe, cmap='jet', vmin = vmin[2], vmax = vmax[2])
            axarr[6][j].axis('off')
            
            
            
            j += 1
        y = 0
        
        vmin[0] = 0
        vmin[1] = 0
        vmax[0] = 1
        vmax[1] = 1
        i = 0
        for c in vmax:
            if i < 2:
                vmax[i] = c
            vmax[i] = vmax[2]
            i +=1
        i = 0
        for c in vmin:
            if i < 2:
                vmin[i] = c
            vmin[i] = vmin[2]
            i +=1

        f.set_size_inches(l*2+4, 12, forward=True)
        cmaps = ['jet']
        for col in range(1):
            for row in range(7):
                ax = axarr[row, col]
                pcm = ax.pcolormesh(np.random.random((0, 0)) * (col + 1),
                                cmap=cmaps[col], vmin = vmin[row], vmax = vmax[row])
                f.colorbar(pcm, ax=ax, location='left')
        for x in nom:  
            axarr[y][0].text(-0.7, 0.45, x,horizontalalignment='right', transform=axarr[y][0].transAxes)
            y += 1

        #plt.savefig('test25_psm.jpg', dpi=600)
        plt.savefig(dossier+'_psm.jpg', dpi=600)
        
        x = np.arange(start, stop+step, step)
        
        
        
        fig = plt.figure(2)
        
        plt.plot(x,y1,'ko-', label = 'raw')
        plt.plot(x,y2,'ro-', label= 'sym')
        plt.ylabel('Puissance relative')
        plt.xlabel('Position relative ($\mu m$)')
        plt.plot(x,y3,'go-', label = 'asym')
        plt.legend()
        #plt.savefig('test25_psm_power1.jpg', dpi=600)
        plt.savefig(dossier+'_psm_power1.jpg', dpi=600)
        fig2 = plt.figure(3)

        
        plt.plot(x,y4,'bo-', label = 'ee')
        plt.plot(x,y5,'yo-', label = 'oo')
        plt.plot(x,y6,'mo-', label = 'eo')
        plt.plot(x,y7,'co-', label = 'oe')
        plt.legend()
        plt.ylabel('Puissance relative')
        plt.xlabel('Position relative ($\mu m$)')
        
        #plt.savefig('test25_psm_power2.jpg', dpi=600)
        plt.savefig(dossier+'_psm_power2.jpg', dpi=600)

        plt.show()

if __name__ == "__main__":
    test = PSM()
    test.psm('test23')
    
    