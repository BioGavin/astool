# Usage: python3 ex_antismash_bgc.py antismash_results antismash_results.xlsx
# antismash_results 文件夹存储着antismash分析的结果文件

from bs4 import BeautifulSoup
import xlsxwriter
import os
import sys

input_dir, output_xlsx = sys.argv[1:3]  # 输入路径
antismash_result_dir = os.path.abspath(input_dir)


workbook = xlsxwriter.Workbook(output_xlsx)
worksheet = workbook.add_worksheet('test')
bold = workbook.add_format({'bold': True})
merge_format = workbook.add_format()
merge_format.set_align('center')
merge_format.set_align('vcenter')
title = ["Name", "Region_Num"]
app = []
worksheet.freeze_panes(1, 0)
org_row = 1
row = 0
col = 0

for i in title:
    worksheet.write(row, col, i, bold)
    col += 1

row = 1

'''
# 该函数不能将A+B和B+A视为同一类
def combineType(data):
    Type = ''
    start = 0
    a = data.select('a')
    for i in a:
        if start == 0:
            Type += i.text
            start += 1
        else:
            Type += '+'+i.text
    return Type
'''


def combineType(data):
    Type = ''
    start = 0
    a = data.select('a')
    if len(a) == 1:
        Type += a[0].text
    else:
        ls = []
        for i in a:
            ls.append(i.text)
        ls.sort()
        Type += "+".join(ls)
    return Type


for file in os.listdir(antismash_result_dir):
    print(file)
    if os.path.exists(os.path.join(antismash_result_dir, file, "index.html")):
        name = file[:15]
        print(name)
        os.chdir(os.path.join(antismash_result_dir, file))
        soup = BeautifulSoup(open('index.html'), features="html.parser")

        try:
            data = soup.select('#overview')[0].select('tbody')[-1]
        except:
            data = 0
        jump = 1
        num = 0
        stat = {}

        if data == 0:
            pass
        else:
            for child in data.children:
                # print(child)
                if len(child) <= 1:
                    continue
                num += 1
                td = child.select('td')
                # print(td)
                Type = combineType(td[1])
                if Type not in app:
                    app.append(Type)
                    worksheet.write(0, len(title) + len(app) - 1, Type, bold)
                if Type not in stat:
                    stat[Type] = 1
                else:
                    stat[Type] += 1

            for i in range(len(app)):
                if app[i] in stat:
                    count = stat[app[i]]
                else:
                    count = 0
                worksheet.write(row, len(title) + i, count)
        worksheet.write(row, 0, name)
        worksheet.write(row, 1, num)
        row += 1
        os.chdir('..')
os.chdir('..')

workbook.close()
