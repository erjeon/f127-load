"""제출 전 문서/이미지 메타데이터 정리. 템플릿 잔재와 도구 흔적을 지운다."""
import sys, os, re, zipfile, shutil, tempfile

CORE_SET = {
 'dc:title':'', 'dc:subject':'', 'dc:creator':'', 'cp:lastModifiedBy':'',
 'dc:description':'', 'cp:keywords':'', 'cp:category':'', 'cp:contentStatus':'',
 'cp:revision':'1',
}
APP_BLANK = ['Company','Manager','Template','HyperlinkBase']

def clean_office(path, title=None):
    tmp = tempfile.mkdtemp()
    with zipfile.ZipFile(path) as z: z.extractall(tmp)
    core = os.path.join(tmp,'docProps','core.xml')
    if os.path.exists(core):
        s = open(core,encoding='utf-8').read()
        for k,v in CORE_SET.items():
            val = title if (k=='dc:title' and title) else v
            s = re.sub(rf'<{k}[^>]*>.*?</{k}>', f'<{k}>{val}</{k}>', s, flags=re.S)
        s = re.sub(r'<cp:lastPrinted>.*?</cp:lastPrinted>','',s,flags=re.S)
        open(core,'w',encoding='utf-8').write(s)
    app = os.path.join(tmp,'docProps','app.xml')
    if os.path.exists(app):
        s = open(app,encoding='utf-8').read()
        for k in APP_BLANK:
            s = re.sub(rf'<{k}>.*?</{k}>', f'<{k}></{k}>', s, flags=re.S)
        s = re.sub(r'<TotalTime>.*?</TotalTime>','<TotalTime>0</TotalTime>',s,flags=re.S)
        open(app,'w',encoding='utf-8').write(s)
    for extra in ('docProps/custom.xml','docProps/thumbnail.jpeg'):
        p=os.path.join(tmp,*extra.split('/'))
        if os.path.exists(p): os.remove(p)
    out = path + '.tmp'
    zf = zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED)
    for root,_,files in os.walk(tmp):
        for f in files:
            fp=os.path.join(root,f)
            zf.write(fp, os.path.relpath(fp,tmp))
    zf.close(); shutil.move(out,path); shutil.rmtree(tmp)
    return True

def clean_png(path):
    from PIL import Image
    im = Image.open(path); data = list(im.getdata()); m = im.mode; sz = im.size
    dpi = im.info.get('dpi',(600,600))
    im2 = Image.new(m,sz); im2.putdata(data)
    im2.save(path, dpi=dpi)          # info 없이 재저장
    return True

if __name__ == '__main__':
    for p in sys.argv[1:]:
        if p.endswith(('.docx','.pptx','.xlsx')):
            clean_office(p, title=None); print("cleaned(office):",os.path.basename(p))
        elif p.lower().endswith('.png'):
            clean_png(p); print("cleaned(png):",os.path.basename(p))
