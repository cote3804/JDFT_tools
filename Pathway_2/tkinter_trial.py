import tkinter as tk

root = tk.Tk()

root.geometry("500x500")
root.title("GUI")

e = tk.Entry(root)
e.pack()

def myClick():
    myLabel = tk.Label(root, text=e.get())
    myLabel.pack()

myButton = tk.Button(root, text="Enter some text dawg", command=myClick)
myButton.pack()

root.mainloop()